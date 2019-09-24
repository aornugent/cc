# Rsq
rsq <- function(output, x, par = "p") {
  
  # Population means
  pop <- extract_long(fit, par, id = "pop") %>%
    spread(par, val) %>%
    rename_at(vars(one_of(par)), ~"mu_pop")
  
  # Plot level effects
  plt <- extract_long(fit, pars = "u", id = "plt") %>%
    spread(par, val)
  
  # Treatment effects
  trt <- extract_long(fit, pars = "beta", id = "ref") %>%
    spread(par, val) %>%
    mutate(trt = as.numeric(ref) + max(x$grp))
  
  left_join(x, pop) %>%
    left_join(plt) %>%
    left_join(trt) %>%
    mutate(beta = if_else(is.na(beta), 1, beta),
           mu_trt = mu_pop * beta,
           mu_plt = mu_pop * u,
           mu_all = mu_trt * u) %>%
    select(sample, abun_std, matches("mu")) %>%
    gather(level, pred, -sample, -abun_std) %>%
    mutate(e = abun_std - pred) %>%
    group_by(level, sample) %>%
    summarise(rsq = var(pred) / (var(pred) + var(e))) %>%
    quantiles("rsq")
}

diagnostics <- function(output) {
  diagnostics <- summary(output$fit)$summary %>%
    as.data.frame() %>%
    rownames_to_column("par") %>%
    select(par, Rhat, n_eff)
  
  message("Worst sampled parameters")
  print(summarise(diagnostics,
            r_hat = max(Rhat),
            n_eff = min(n_eff)))
  
  posterior_summary(diagnostics)
}

extract_long <- function(fit, pars, id = NA) {
  divergent <- get_sampler_params(fit, inc_warmup = F) %>%
    map_df(., as.data.frame) %>%
    rename_all(~ gsub("_", "", .)) %>%
    select(divergent) %>%
    mutate_all(~ if_else(. == 1, T, F)) %>%
    rowid_to_column("sample")

  samples <- extract(fit, pars = pars, 
                     permuted = F,
                     inc_warmup = F) %>%
    as.data.frame() %>%
    rowid_to_column("iteration") %>%
    gather(chain, val, -iteration) %>% 
    separate(chain, c("chain", "par"), sep = "\\.") %>%
    spread(par, val) %>%
    rowid_to_column("sample") %>%
    gather(par, val, -chain, -iteration, -sample) %>%
    left_join(divergent)
  
  if(!is.na(id)) {
    samples <- separate(samples, par, c("par", id),
                        sep = "\\[|\\]", extra = "drop") %>%
    mutate_at(vars(matches(id)), as.numeric)
  }

  samples <- spread(samples, par, val)
  return(samples)
}

quantiles <- function(df, par) {
  summarise_at(df, vars(one_of(par)),
              .funs = list(mean = ~ mean(.), 
                  low = ~ quantile(., 0.025),
                  high = ~ quantile(., 0.975)))
}

posterior <- function(output) {
  
  x <- get_x(get_y(output$fields, quietly = T))
  
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
  
  pop <- extract_long(output$fit, pars = pars, id =  "pop") %>%
    left_join(get_p(x)) 
  
  # Plot level effects
  plt <- extract_long(output$fit, pars = "u", id = "plt") %>%
    left_join(get_u(x))
  
  # Treatment effects
  trt <- extract_long(output$fit, pars = "beta", id = "ref") %>%
    mutate(trt = as.numeric(ref) + max(x$grp)) %>%
    left_join(get_g(x))
  
  # Time steps
  ts <- output$data_list$T
  
  post <- left_join(pop, plt) %>%
      left_join(trt) %>%
      left_join(x)
      
  rm(list = c("output", "pop", "plt", "trt"))
  gc()
  
  if (grepl("f2", model)) {
    growth <- extract_long(output$fit, pars = "rK", id = "grp") %>%
      left_join(get_p(x))
    
    post = left_join(post, growth)
    
  } else if (grepl("f3", model)) {
    pred = mutate(pop, pred = beta_curve(p0, pK, tK, tmax, t, n()))
  }
  
  return(pred)
}
