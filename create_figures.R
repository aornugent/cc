theme_set(
  theme_bw() + 
    theme(aspect.ratio = 1,
          text = element_text(size = 12))
)

posterior_summary <- function(diagnostics) {
  gather(diagnostics, meas, val, -par) %>%
    ggplot(., aes(x = val)) +
    geom_histogram(fill = "blue", 
                   alpha = 0.2) +
    coord_cartesian(expand = F) +
    facet_grid(~ meas, scales = "free")
}

figure_one <- function(save = T) {
  x <- data.frame(year = c(1993, 1994, 1999, 2000, 2017),
                  intercept = c(3.5, 3.5, 2.5, 2.5, 3.5))
  
  y <- get_y(92) %>%
    filter(group != "Other") %>%
    left_join(x)
  
  
  p <- ggplot(y, aes(x = rest, y = abun_std)) +
    geom_boxplot(aes(fill = group), width = 0.8) + 
    geom_vline(aes(xintercept = intercept), linetype = "dashed") +
    guides(fill = FALSE) + 
    facet_grid(group ~ year, scales = "free") +
    labs(x = "", y = "Abundance (standardised)")+
    theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.3)) 
  
  if(save)
    ggsave(p, filename = "figs/figure_one.png", device = "png",
           dpi = 600, width = 22, height = 18, units = "cm")
  
  return(p)
}

figure_two <- function(x, output, stat = "mode", 
                       grouped_rK = F, save = T) {
  message("Loading model")
  model <- output$model
  
  expose_stan_functions(output$model_file)
  
  pars <- switch(model,
                 f1 = c("p"),
                 f1b = c("p"),
                 f2 = c("p0", "pK"),
                 f2b = c("p0", "pK"),
                 f3 = c("p0", "pK", "tK", "tmax"),
                 f3b = c("p0", "pK", "tK", "tmax"))
  
  time <- get_t(x)
  
  pop <- extract_long(output$fit, pars = pars, id =  "pop") %>%
    left_join(get_p(x)) %>%
    group_by(sample, pop) %>%
    mutate(t = list(0:max(time$t))) %>%
    unnest() 
  
  if(grepl("f1", model)) {
    pop = rename_at(pop, vars(pars), ~ "pred")
  } else if (grepl("f2", model)) {
    if(grouped_rK)
      growth <- extract_long(output$fit, pars = "rK", id = "grp") %>%
        left_join(get_g(x))
    else {
      growth <- extract_long(output$fit, pars = "rK", id = "pop") %>%
        left_join(get_p(x))
    }
    
    sigma <- extract_long(output$fit, pars = "sigma_e", id = "grp") %>%
      left_join(get_p(x))
    
    pop = left_join(pop, growth) %>%
      left_join(sigma) %>%
      mutate(abun = gompertz_curve(p0, pK, rK, t, n()),
             mode = abun * (1 + sigma_e^2 / abun^2)^-1.5)
  
    } else if (grepl("f3", model)) {
    pop = mutate(pop, pred = beta_curve(p0, pK, tK, tmax, t, n()))
  }
  
  pred <- group_by(pop, pop, t) %>%
    quantiles(stat) %>%
    left_join(get_p(x)) %>%
    mutate(year = t + ab,
           field = factor(field, labels = c("F1", "F2", "F0"),
                          levels = c(44, 53, 92))) %>%
    filter(year %in% time$year) %>%
    filter(group != "Other")
  
  x <- mutate(x, hide = if_else(rest %in% c("R-", "C-"), F, T),
              field = factor(field, labels = c("F1", "F2", "F0"),
                             levels = c(44, 53, 92))) %>%
    filter(group != "Other")
  
  delta <- extract_long(output$fit, pars = "delta", id =  "grp") %>%
    left_join(get_g(x)) %>%
    group_by(group) %>%
    summarise(label = sprintf("~ delta == %0.2f", mean(delta))) %>%
    mutate(x = 2010,
           y = 5.3,
           field = "F0") %>%
    filter(group != "Other")
  
  
  p <- ggplot(pred, aes(x = year, y = mean)) +
    geom_line(data = x, size = 0.7, 
              aes(y = abun_std, group = plot_id, alpha = hide)) +
    geom_ribbon(aes(ymin = low, ymax = high),
                alpha = 0.7, fill = "royalblue2") +
    geom_line(aes(y = mean), color = "white",
              size = 1.2, linetype = "dashed") +
    geom_text(data = delta, aes(x = x, y = y, label = label), parse = T) +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    coord_cartesian(expand = F) + 
    facet_grid(field ~ group, scales = "free") +
    labs(x = "", y = "Abundance (standardised)", 
         fill = "") +
    guides(linetype = F, alpha = F) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.3),
          legend.position = "bottom") 
  
  if(save)
    ggsave(p, filename = "figs/figure_two.png",
           dpi = 600, width = 18, units = "cm")
  
  return(p)
}
  

#Treatment effects
figure_two_b <- function(x, output, save = T) {
  beta <- as.data.frame(output$fit, pars = "beta") %>%
    gather(par, val) %>%
    separate(par, c("par", "ref"),
             sep = "\\[|\\]", extra = "drop") %>%
    mutate(trt = as.numeric(ref) + max(x$grp)) %>%
    group_by(par, ref, trt) %>%
    summarise(mean = mean(val),
              low = quantile(val, 0.025),
              high = quantile(val, 0.975)) %>%
    left_join(get_tr(x)) %>%
    filter(group != "Other")


  p <- ggplot(beta, aes(x = rest, y = mean, ymin = low, ymax = high)) +
    geom_point() +
    geom_errorbar(width = 0) +
    geom_hline(aes(yintercept = 1)) +
    facet_grid(~ group) +
    labs(x = "", y = "Effect size") + 
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.3),
          legend.position = "bottom") 
  
  if(save)
    ggsave(p, filename = "figs/figure_two_b.png",
           dpi = 600, width = 16, height = 6, units = "cm")
  
  return(p)
}

figure_three <- function(x, output, stat = "abun", level = "grp",
                         save = T) {
  expose_stan_functions(output$model_file)
  time <- get_t(x)
  
  if(level == "grp") {
    pars = c("mu_p0", "mu_pK", "mu_rK", "sigma_e")
    grp <- extract_long(output$fit, pars, id =  "grp") %>%
      left_join(get_g(x)) %>%
      group_by(sample, group) %>%
      mutate(t = list(0:max(time$t))) %>%
      unnest() %>%
      mutate(abun = gompertz_curve(mu_p0, mu_pK, mu_rK, t, n()),
             mode = abun * (1 + sigma_e^2 / abun^2)^-1.5)
  } else if(level == "pop") {
    pars = c("p0", "pK", "rK")
    
    sigma <- extract_long(output$fit, pars = "sigma_e", id = "grp") %>%
      left_join(get_p(x))
    
    grp <- extract_long(output$fit, pars, id =  "pop") %>%
      left_join(get_p(x)) %>%
      left_join(sigma) %>%
      group_by(sample, pop) %>%
      mutate(t = list(0:max(time$t))) %>%
      unnest() %>%
      mutate(abun = gompertz_curve(p0, pK, rK, t, n()),
             mode = abun * (1 + sigma_e^2 / abun^2)^-1.5)
  }
    
  pred <- group_by(grp, group, t) %>%
    quantiles(stat) %>%
    left_join(get_p(x)) %>%
    filter(group != "Other")
  
  v <- left_join(x, time) %>%
    mutate(hide = if_else(rest %in% c("R-", "C-"), F, T)) %>%
    filter(group != "Other")
   
  delta <- extract_long(output$fit, pars = "delta", id =  "grp") %>%
    left_join(get_g(x)) %>%
    group_by(group) %>%
    summarise(label = sprintf("~ delta == %0.2f", mean(delta))) %>%
    mutate(x = 70,
           y = c(8.2, 9.5, 7.7, 8.5, 11.3),
           field = 44) %>%
    filter(group != "Other")
  
  p <- ggplot(v, aes(x = t)) +
    geom_line(aes(y = abun_std, group = plot_id, colour = ab, alpha = hide)) +
    geom_ribbon(data = pred, aes(ymin = low, ymax = high), 
                fill = "goldenrod2", alpha = 0.6) +
    geom_line(data = filter(v, field == 92), color = "white",
              aes(y = abun_std, group = plot_id, alpha = hide)) + 
    geom_line(data = pred, aes(y = mean), color = "yellow",
              size = .8, linetype = "dashed") +
    geom_text(data = delta, aes(x = x, y = y, label = label), parse = T) +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    coord_cartesian(expand = F) +
    facet_grid(group ~ ., scales = "free") +
    labs(x = "Years since abandonment",
         y = "Abundance (standardised)", 
         fill = "", colour = "") +
    guides(linetype = F, alpha = F) +
    theme_bw() +
    theme(aspect.ratio = 0.7,
      legend.position = "right") 
  
  if(save)
    ggsave(p, filename = "figs/figure_three.png",
           dpi = 600, width = 10, height = 17, units = "cm")
  
  return(p)
}

figure_four <- function(x, output) {
  
  
  x <- mutate(x, hide = if_else(rest %in% c("R-", "C-"), F, T)) %>%
    left_join(time)
  
  ggplot(pop, aes(x = t, y = mean)) +
    geom_line(data = x, size = 0.7, 
              aes(y = abun_std, group = plot_id, 
                  alpha = hide, color = ab)) +
    geom_ribbon(aes(ymin = low, ymax = high, group = field), 
                fill = "white", color = "black", alpha = 0.4) +
    geom_line(aes(y = mean, group = field),
              color = "red", size = .8, linetype = "dashed") +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    labs(x = "", y = "Abundance (standardised)", 
         fill = "") +
    guides(linetype = F, alpha = F) +
    facet_grid(~ group) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.3),
          legend.position = "bottom")
}

# Variety of lognormals, moments and transformations
figure_zero_log <- function() {
  sim <- data.frame(group = LETTERS[1:2],
                     log_mu = 1, 
                     log_sigma = c(0.6, 2)) %>%
    mutate(mean = exp(log_mu + 0.5 * log_sigma^2),
           mode = exp(log_mu - log_sigma^2)) %>%
    mutate(x = list(seq(0.1, 20, length.out = 100)),
           log_x = list(seq(-1, 4, length.out = 100))) %>%
    unnest() %>%
    mutate(dens = dlnorm(x, log_mu, log_sigma),
           log_dens = dnorm(log_x, log_mu, log_sigma)) %>%
    filter(dens > 1e-4)
  
  label <- list(bquote(sigma == 0.6),
                bquote(sigma == 2.0))
  
  ggplot(sim, aes(x = x, y = dens, group = group,
                  colour = group)) +
    geom_line(size = 1.1) +
    scale_colour_manual(values = c("blue4", "royalblue2")) +
    guides(color = F)

  ggplot(sim, aes(x = log_x, y = log_dens, group = group,
                  colour = group)) +
    geom_line(size = 1.1) +
    scale_colour_manual(values = c("blue4", "royalblue2"),
                        label = label) +
    labs(color = "",
         y = "Density", x = "log(x)") +
    coord_flip() +
    scale_x_reverse() +
    theme(legend.position = c(0.9, 0.1))
  
  ggplot(sim, aes(x = exp(log_x), y = log_x, group = group)) +
    geom_line(size = 1.1, color = "black") +
    scale_y_reverse() +
    guides(color = F) + 
    labs(color = "", y = "", x = "x")
}

figure_zero <- function(output, x,
                        pars = c("mu_p0", "mu_pK", "mu_rK", "sigma_e"),
                        n = 80, save = T) {
  
  expose_stan_functions(output$model_file)
  
  ts = output$data_list$T
  
  stat <- extract_long(output$fit, pars, "grp") %>%
    left_join(get_g(x)) %>%
    group_by(group) %>%
    quantiles(pars) %>%
    gather(quantile, val, -group) %>%
    separate(quantile, c("par", "quantile"), sep = "_(?=[^_]+$)") %>%
    spread(par, val) %>%
    mutate_at(vars(quantile), ~fct_relevel(., "low", "mean", "high")) %>%
    group_by(group, quantile) %>%
    mutate(t = list(seq(1, ts, len = n))) %>%
    unnest() %>%
    mutate(mu_p = gompertz_curve(mu_p0, mu_pK, mu_rK, t, n())) %>%
    mutate(log_sigma = sqrt(log(1 + sigma_e^2 / mu_p^2)),
           log_mu = log(mu_p) - 0.5 * log_sigma^2)
  
  sim <- select(stat, group, quantile, t, log_mu, log_sigma) %>%
    mutate(x = list(seq(0.1, 10, length.out = 1000))) %>%
    unnest() %>%
    mutate(dens = dlnorm(x, log_mu, log_sigma))
  
  p <- ggplot(sim, aes(x = x)) +
    geom_line(aes(y = dens, group = t, color = t),
              alpha = 0.1, size = 1.1) +
    facet_grid(quantile ~ group, scales = "free") +
    stat_density(data = filter(x, abun_std >= 0, rest == "R-"), 
                 aes(x = abun_std)) +
    coord_cartesian(ylim = c(0, 2.0), expand = F) +
    scale_colour_gradient(low = "darkred", high = "red") +
    labs(x = "Abundance (standardised)",
         y = "Density") +
    theme(legend.position = "bottom")
  
  if(save)
    ggsave(p, filename = "figs/figure_appendix_one.png",
           dpi = 600, width = 18, units = "cm")
  
  return(p)
  
  return(p)
}



posterior_density <- function(output, pars, id = "id") {
  extract_long(output$fit, pars, id) %>%
    gather(par, val, 
           -c("sample", "iteration", "chain", "id", "divergent")) %>%
    ggplot(., aes(x = val, fill = factor(divergent))) +
    geom_histogram(position = "identity") +
    scale_fill_manual(values = c("darkgrey", "darkred")) +
    labs(fill = "Divergent samples", x = "") +
    facet_wrap(~ par, scales = "free")
}