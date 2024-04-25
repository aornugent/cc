library(tidyverse)
library(rstan)

source("cc/run_model.R")
source("cc/check_model.R")

# Lajeunesse 2011, Eq. 3
agg_rr <- function(E, V) {
  X = rep(1, nrow(V))
  V_inv = solve(V)
  
  RR_var = (X %*% V_inv %*% X)^-1
  RR = RR_var * (X %*% V_inv %*% E)
  
  return(RR_var)
}

rr <- function(df) {
  control_response <- filter(df, rest == "R-")

  # Lajeunesse 2011, Eq. 1
  response_ratios <- filter(df, rest != "R-" & rest != "B-") %>%
    left_join(control_response, by = "group", suffix = c("_t", "_c")) %>%
    mutate(rr = log(mean_t / mean_c))
  
  # Lajeunesse 2011, Eq. 1, 8, 3
  rr <- group_by(response_ratios, group) %>%
    summarise(E = list(rr),
              V = list(diag(sd_t^2 / (n_t * mean_t^2)) + 
                         sd_c^2 / (n_c * mean_c^2)),
              RR_var = agg_rr(E, V)) %>%
    left_join(select(response_ratios, group, n_c, rest_t, rr)) %>%
    spread(rest_t, rr) %>%
    select(-E, -V) %>%
    mutate(RR_sd = sqrt(RR_var))
  
  return(rr)
}

# Local scale
x <- get_y(fields = 92) %>%
  get_x(.)

group_response <- filter(x, year == 2017) %>%
  group_by(group, rest) %>%
  summarise(n = n(),
            mean = mean(abun_std),
            sd = sd(abun_std))

rr_local <- rr(group_response)


# Spatial scale
post <- function(fit, x, year = 52) {
  post_tr <- extract_long(fit, pars = "beta", id = "ref") %>%
    mutate(trt = as.numeric(ref) + max(x$grp)) %>%
    left_join(get_tr(x))
  
  # post_plt <- extract_long(fit, pars = "u", id = "plt") %>%
  #   left_join(get_u(x))
  
  post_mean <- extract_long(fit, pars = c("p0", "pK", "rK"), id = "pop") %>%
    left_join(get_p(x)) %>%
    left_join(get_tr(x)) %>%
    # left_join(get_u(x)) %>%
    left_join(post_tr) %>%
    # left_join(post_plt) %>%
    mutate(beta = if_else(is.na(beta), 1, beta)) %>%
    mutate(x = gompertz_curve(p0, pK * beta, rK, rep(year, n()), n()))
  
  return(post_mean)
}


load("cc/f2/field_mean2f3v1censored.fit.Rdata")
expose_stan_functions(output$model_file) 

x <- get_y(fields = output$fields) %>%
  get_x(.)

post_spatial <- post(output$fit, x)

group_response <- select(x, field, rest, group, plot_id) %>%
  unique() %>%
  left_join(post_spatial) %>%
  group_by(group, rest) %>%
  summarise(mean = mean(x),
            sd = sd(x),
            n = n_distinct(plot_id),
            n_samples = n())

rr_spatial <- rr(group_response)

# Temporal scale
load("cc/f2/field_mean2f22v1censored.fit.Rdata")
x <- get_y(output$fields) %>%
  get_x(.)

post_temporal <- post(output$fit, x)
group_response <- select(x, field, rest, group, plot_id) %>%
  unique() %>%
  left_join(post_temporal) %>%
  group_by(group, rest) %>%
  summarise(mean = mean(x),
            sd = sd(x),
            n = n_distinct(plot_id),
            n_samples = n())

rr_temporal <- rr(group_response)


rr <- bind_rows(rr_local, rr_spatial) %>%
  bind_rows(rr_temporal)

write_csv(x = rr, path = "cc/response_ratios2.csv")
