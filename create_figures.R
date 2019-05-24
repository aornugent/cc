theme_set(
  theme_bw() + 
    theme(aspect.ratio = 1)
)

posterior_summary <- function(diagnostics) {
  gather(diagnostics, meas, val, -par) %>%
    ggplot(., aes(x = val)) +
    geom_histogram(fill = "blue", 
                   alpha = 0.2) +
    coord_cartesian(expand = F) +
    facet_grid(~ meas, scales = "free")
}

figure_one <- function(y) {
  # Figure 1
  ggplot(y, aes(x = rest, y = abun_std, 
                group = rest, fill = group)) +
    geom_boxplot(width = 0.8) + 
    facet_grid(group ~ year, scales = "free_x") +
    labs(x = "", y = "Abundance (standardised)", 
         fill = "", title = "Restoration experiment") +
    theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.3),
          legend.position = "bottom") 
}

figure_two <- function(x, output, stat = "mode") {
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
    growth <- extract_long(output$fit, pars = "rK", id = "pop") %>%
      left_join(get_p(x))
    
    sigma <- extract_long(output$fit, pars = "sigma_e_sq", id = "grp") %>%
      left_join(get_p(x))
    
    pop = left_join(pop, growth) %>%
      left_join(sigma) %>%
      mutate(mean = gompertz_curve(p0, pK, rK, t, n()),
             mode = mean / (sigma_e_sq / mean^2 + 1)^1.5)
  
    } else if (grepl("f3", model)) {
    pop = mutate(pop, pred = beta_curve(p0, pK, tK, tmax, t, n()))
  }
  
  pred <- group_by(pop, pop, t) %>%
    quantiles(stat) %>%
    left_join(get_p(x)) %>%
    mutate(year = t + ab) %>%
    filter(year %in% time$year)
  
  x <- mutate(x, hide = if_else(rest %in% c("R-", "C-"), F, T))
  
  ggplot(pred, aes(x = year, y = mean)) +
    geom_line(data = x, size = 0.7, 
              aes(y = abun_std, group = plot_id, alpha = hide)) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = group), 
                color = "black", alpha = 0.4) +
    geom_line(aes(y = mean), color = "red",
              size = .8, linetype = "dashed") +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    facet_grid(group ~ field, scales = "free") +
    labs(x = "", y = "Abundance (standardised)", 
         fill = "") +
    guides(linetype = F, alpha = F) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.3),
          legend.position = "bottom") 
}
  

#Treatment effects
figure_two_b <- function(x, output) {
  beta <- as.data.frame(output$fit, pars = "beta") %>%
    gather(par, val) %>%
    separate(par, c("par", "ref"),
             sep = "\\[|\\]", extra = "drop") %>%
    mutate(trt = as.numeric(ref) + max(x$grp)) %>%
    group_by(par, ref, trt) %>%
    summarise(mean = mean(val),
              low = quantile(val, 0.025),
              high = quantile(val, 0.975)) %>%
    left_join(get_g(x))


  ggplot(beta, aes(x = rest, y = mean, ymin = low, ymax = high)) +
    geom_point() +
    geom_errorbar(width = 0) +
    geom_hline(aes(yintercept = 1)) +
    facet_grid(group ~ .) +
    labs(x = "", y = "Effect size (%)") + 
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.3),
          legend.position = "bottom") 
}

figure_three <- function(x, output, par = "p", all = F) {
  sigma <- extract_long(output$fit, pars = "sigma_e_sq", id = "grp") %>%
    spread(par, val) 
  
  pop <- extract_long(output$fit, par, id =  "pop") %>%
    spread(par, val) %>%
    rename_at(vars(one_of(par)), ~ "mu_pop") %>%
    left_join(get_p(x)) %>%
    left_join(sigma) %>%
    mutate(log_sigma_sq = log(1 + sigma_e_sq / mu_pop^2),
           log_mu = log(mu_pop) - 0.5 * log_sigma_sq)
  
  dens <- group_by(pop, field, group) %>%
    summarise_at(vars("log_mu", "log_sigma_sq"), mean) %>%
    group_by(field, group) %>%
    mutate(abun_std = list(seq(0, 6, len = 1000))) %>%
    unnest() %>%
    mutate(dens = dlnorm(abun_std, log_mu, sqrt(log_sigma_sq)))
  
  pop_mean <- group_by(pop, field, group) %>%
    summarise(mu_pop = mean(mu_pop))
  
  filter(x, grepl(if_else(all, ".", "R-|C-"), rest)) %>%
    ggplot(., aes(x = abun_std, fill = group)) +
      geom_density(alpha = 0.2) + 
      geom_line(data = dens, aes(y = dens), color = "red") + 
      geom_vline(data = pop_mean, aes(xintercept = mu_pop)) +
      coord_cartesian(ylim = c(0, 2), xlim = c(0, 5),
                      expand = F) +
      facet_grid(field ~ group) +
      labs(x = "Abundance (standardised)", y = "Density",
           fill = "") +
      guides(linetype = F, alpha = F) +
      theme_bw() +
      theme(aspect.ratio = 1,
            legend.position = "bottom") 
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
  stat <- data.frame(log_mu = 1, 
                     log_sigma = seq(0.4, 0.8, len = 100)) %>%
    mutate(mean = exp(log_mu + 0.5 * log_sigma^2),
           mode = exp(log_mu - log_sigma^2),
           median = exp(log_mu)) %>%
    gather(stat, x, mean:median) %>%
    mutate(dens = dlnorm(x, log_mu, log_sigma))
  
  sim <- select(stat, log_mu, log_sigma) %>%
    unique() %>%
    filter(log_sigma %in% c(min(log_sigma),  max(log_sigma))) %>%
    mutate(x = list(seq(0.1, 10, length.out = 1000))) %>%
    unnest() %>%
    mutate(dens = dlnorm(x, log_mu, log_sigma))
  
  ggplot(sim, aes(x = x, y = dens, group = log_sigma,
                  alpha = log_sigma)) +
    geom_line(color = "blue", size = 1.1) +
    geom_segment(data = stat, aes(xend = x, yend = 0,
                                  color = stat), size = 1.1) +
    scale_alpha_continuous(range = c(0.1, 0.8)) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  log_sim <- mutate(sim, log_x = log(x),
                    log_dens = dnorm(log_x, log_mu, log_sigma))
  
  ggplot(log_sim, aes(x = x, y = log_x, group = log_sigma)) +
    geom_line() +
    scale_y_reverse() +
    scale_alpha_continuous(range = c(0.1, 0.8)) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  ggplot(log_sim, aes(x = log_x, y = log_dens,
                      group = log_sigma, alpha = log_sigma)) +
    geom_line(color = "blue", size = 1.1) +
    geom_line() +
    scale_x_reverse() +
    coord_flip() +
    scale_alpha_continuous(range = c(0.1, 0.8)) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
}

figure_zero <- function(fit, pars = c("mu_pK", "sigma_e")) {
  stat <- extract_long(fit, pars, "grp") %>%
    group_by(grp) %>%
    quantiles(pars) %>%
    gather(quantile, val, -grp) %>%
    separate(quantile, c("par", "quantile"), sep = "_(?=[^_]+$)") %>%
    spread(par, val) %>%
    mutate_at(vars(quantile), ~ fct_relevel(quantile, "mean")) %>%
    rename_at(vars(pars), ~ c("mu", "sigma")) %>%
    mutate(log_sigma = sqrt(log(1 + sigma^2 / mu^2)),
           log_mu = log(mu) - 0.5 * log_sigma^2,
           mean = exp(log_mu + log_sigma^2 / 2),
           mode = exp(log_mu - log_sigma^2),
           median = exp(log_mu)) %>%
    gather(stat, x, mean:median) %>%
    mutate(dens = dlnorm(x, log_mu, log_sigma))
  
  sim <- select(stat, grp, quantile, log_mu, log_sigma) %>%
    unique() %>%
    mutate(x = list(seq(0.2, 5, length.out = 1000))) %>%
    unnest() %>%
    mutate(dens = dlnorm(x, log_mu, log_sigma))
  
  ggplot(sim, aes(x = x, y = dens, group = paste(grp, quantile),
                  color = factor(grp), linetype = quantile)) +
    geom_line(size = 1.1) +
    #geom_vline(data = stat, aes(xintercept = mu_pK)) +
    #geom_segment(data = stat, aes(xend = x, yend = 0,
    #                              color = stat), size = 1.1) +
    facet_grid(~ grp) +
    scale_alpha_continuous(range = c(0.1, 0.8)) +
    theme_bw() +
    theme(aspect.ratio = 1)
}