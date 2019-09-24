#' utils

run_models <- function(model = "f1", fields = 92, 
                       censored = T, variant = 1,
                       chains = 4, iter = 400, init_r = 0.5) {

  output <- list(
    model = model,
    fields = fields,
    censored = censored,
    model_file = switch(model,
      f1 = "cc/f1/field_mean1.stan",
      f1b = "cc/f1/field_mean1b.stan",      
      f2 = "cc/f2/field_mean2.stan",
      f2b = "cc/f2/field_mean2b.stan",
      f3 = "cc/f3/field_mean3.stan"),
    data_list = format_data(fields)
  )
  
  if(model != "f3")
    variant = 1
  
  output$data_list$model_variant = variant
  output$data_list$cen = if_else(censored, 1, 0)

  output$fit <- stan(file = output$model_file,
                     data = output$data_list,
                     chains = chains,
                     iter = iter,
                     init_r = init_r)
  
  censored <- switch(censored,
                     T = "censored",
                     F = "uncensored")
  
  n_fields <- output$data_list$N_pop / output$data_list$N_grp 
        

  output$ouput_file <- paste0(gsub("\\.stan", "", output$model_file), 
                              "f", n_fields, "v", variant, censored, ".fit.Rdata")
  
  save(file = output$ouput_file, "output")
  message("Model saved")
  
  str(output, max.level = 1)
  return(output)
}


format_data <- function(fields) {

  y <- get_y(fields, quietly = T)
  L <- get_L(y)
  x <- get_x(y)
  t <- get_t(x)
  
  data_list <- list(
    T = max(t$t),
    N = nrow(x),
    N_grp = max(x$grp),
    N_trt = max(x$trt),
    N_plt = max(x$plt),
    N_pop = max(x$pop),
    N_meas = max(x$meas),
    N_mis = sum(x$mis),
    t = t$t,
    y = x$abun_std,
    gm = t$gm,
    L = L$min,
    m1 = t$m[t$m == t$m1],
    m = t$m[t$m != t$m1],
    m_m1 = t$mm1[t$m != t$m1],
    m_obs = t$m[t$mis == 0],
    m_mis = t$m[t$mis == 1],
    grp = x$grp,
    grp_plt = nested(x$plt, x$grp),
    grp_pop = nested(x$pop, x$grp),
    grp_trt = nested(x$trt[x$trt > max(x$grp)], 
                     x$grp[x$trt > max(x$grp)]),
    trt = x$trt,
    plt = x$plt,
    pop = x$pop,
    meas = x$meas
  )
  
  return(data_list)
}

load_data <- function(fields) {
  # CC chronosequence + restoration experiment.
  # Chronosequence plots have burned and unburned treatments.
  # Restoration has 7 different contrasts, including negative controls.
  # Consider C- same as R-, ignoring quadrat level differences.
  plots <- read_csv("data/CC_plot_data.csv", guess_max = 1e4) %>%
    mutate(rest = case_when(burned == 0 ~ "R-",
                            burned == 1 ~ "RB-",
                            T ~ restoration93),
           rest = factor(rest,
                         labels = c("C-", "B-", "C-", "C+",
                                     "B+", "HB+", "BR+", "BRN+", "HBR+"), 
                         levels = c("R-", "RB-", "C-", "C+",
                                      "B", "HB", "BR", "BRN", "HBR")),
           treat = factor(treatment99, levels = c("C", "N", "B", "N/B"))) %>%
    select(plot_id, field, block, rest, yearabandoned, quad_id, treat) %>%
    distinct()
  
  # Select three fields of similar abandonment age. 
  # Bray-Curtis says 42 & 53 most similar of all CC fields to C- plots in 1993. 
  # Group by origin and functional group. 
  # Note: Poa pratensis considered 'Introduced' in this study as origin 
  #       unknown but a weed in Minnesota. 
  dat <- read_csv("data/CC_species_data.csv", 
                  col_types = "ciiccdddddcccccccc") %>%
    select(-count, -pres) %>%
    gather(meas, val, biomass:pcover) %>%
    filter(val > 0) %>%
    left_join(plots) %>%
    filter(field %in% fields) %>%
    filter(!grepl("ground|litter", species_id)) %>%
    mutate(origin = if_else(sp == "Poa_pra", "Introduced", origin)) %>%
    mutate(group = case_when(lifeform == "Grass" & 
                               origin == "Introduced" ~ "Non-native grasses",
                             lifeform == "Grass" &
                               origin == "Native" ~ "Native grasses",
                             func_grp == "LF" &
                               origin == "Introduced" ~ "Non-native forbs",
                             func_grp == "F" &
                               origin == "Introduced" ~ "Non-native forbs",
                             func_grp == "LF" & origin == "Native" ~ "Native forbs",
                             func_grp == "F" & origin == "Native" ~ "Native forbs",
                             func_grp == "W" ~ "Other",
                             T ~ "Other")) %>%
    mutate(group = factor(group, 
                          levels = c("Non-native grasses", 
                                     "Native grasses",
                                     "Non-native forbs",
                                     "Native forbs", 
                                     "Other")))
  
}

get_y <- function(fields, quietly = F) {
  withCallingHandlers(
    dat <- load_data(fields),
    message = function(message){if(quietly) invokeRestart("muffleMessage")}
  )
  
  # Aggregte groups at plot level in each year.
  # Measurement techniques vary significantly between surveys,
  # Standardise by meaurement type and year to maintain relative abundance 
  # without converting to proportions. 
  # Add zeros to plots where no abundance observed.
  
  message("Generating data")
  y <- group_by(dat, field, year, group) %>%
    group_by(field, year, plot_id, meas, rest, ab = yearabandoned, group) %>%
    summarise(abun = sum(val)) %>%
    complete(group, fill = list(abun = 0)) %>%
    group_by(meas, year) %>%
    mutate(sd = sd(abun),
           abun_std = abun / sd) %>%
    ungroup()
}

get_L <- function(y) {
  # Find detection limits of each survey
  L <- group_by(y, year) %>%
    summarise(min = min(abun_std[abun_std > 0])) %>%
    ungroup() %>%
    mutate(meas = as.numeric(factor(paste(year)))) %>%
    arrange(year)
}

get_x <- function(y) {
  # Generate unique IDs for each group in each plot, treatment.
  # Reference classes rely on R-/C- being first factor level.
  # Generate ID for each survey, indicator for non-detections.
  x <- y %>% 
    mutate_at(vars(rest), ~ fct_recode(., `R-` = "C-")) %>%
    mutate(grp = as.numeric(factor(group)),
           plt = as.numeric(factor(paste(as.numeric(group), plot_id))),
           pop = as.numeric(factor(paste(as.numeric(group), field))),
           trt = grp + (as.numeric(rest) - 1) * max(grp),
           meas = as.numeric(factor(year)),
           mis = if_else(abun_std == 0, 1, 0)) %>%
    arrange(plt, year)
}

get_t <- function(x) {
  # Calculate time since abandonment, time between surveys for each plot.
  # Generate indices for first, lagged measurements
  t = select(x, ab, year, plt, grp, mis) %>%
    mutate(t = year - ab,
           m = row_number()) %>%
    group_by(plt) %>%
    mutate(m1 = min(m),
           mm1 = lag(m, default = 0),
           gm = year - lag(year, default = unique(ab))) %>%
    ungroup() %>%
    arrange(plt, year)
  
}

nested <- function(x, y){
  nestbl <- table(factor(x), factor(y))
  idx <- unname(apply(nestbl, 1, function(x) which(x > 0)))
  return(idx)
}

get_p <- function(x) {
  select(x, field, group, pop, grp, ab) %>%
    unique()
}

get_u <- function(x) {
  select(x, field, plot_id, group, plt) %>%
    unique()
}

get_g <- function(x) {
  select(x, grp, group) %>%
    unique()
}

get_tr <- function(x) {
  select(x, trt, rest, group) %>%
    unique()
}