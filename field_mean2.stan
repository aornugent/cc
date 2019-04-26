data {
  int<lower=1> T;       // timespan
  int<lower=1> N;       // total obs
  int<lower=1> N_grp;   // groups
  int<lower=1> N_trt;   // contrasts
  int<lower=1> N_plt;   // experiments
  int<lower=1> N_pop;   // populations
  int<lower=1> N_obs;   // observers
  
  vector[N] y_obs;      // abundance
  vector[N] g_m;        // time between obs.
  
  vector[N_obs] obs_sd; // meas. scale
  int<lower=1, upper=N_obs> obs[N]; // meas. index
  
  // obs indices
  int<lower=1, upper=N> t_1[N_plt];
  int<lower=1, upper=N> t[N - N_plt];
  int<lower=1, upper=N> t_m1[N - N_plt];

  // plot, pop., group & treatment indices
  int<lower=1, upper=N_grp> grp[N];
  int<lower=1, upper=N_pop> pop[N]; 
  int<lower=1, upper=N_plt> plt[N];
  int<lower=1, upper=N_grp> grp_p[N_plt];
  int<lower=1, upper=N_trt> trt[N];
  int<lower=1, upper=N_grp> grp_t[N_trt - N_grp];
}
transformed data{
  // Standardise and offset for log transformation
  vector[N] log_y_std = log((y_obs ./ obs_sd[obs]) + obs_sd[obs]);
}
parameters{
  // intercepts
  vector[N_pop] log_p0_raw;
  
  // endpoints
  vector[N_pop] log_p_raw;
  real<lower=0> sigma_log_p;

  // predictors ordered so R- x N_grp at start
  vector[N_trt - N_grp] log_beta;
  vector[N_grp] mu_log_b;
  vector<lower=0>[N_grp] sigma_log_b;
  
  // correlation
  real<lower=0, upper=1> delta;
  
  // random effects
  vector[N_plt] log_u_raw;
  vector<lower=0>[N_grp] sigma_log_u;
  
  // obs. error
  vector<lower=0>[N_grp] sigma_log_e;
}
model{
  {
    vector[N_pop] log_p0;
    vector[N_pop] log_p;
    vector[N_plt] log_u;
    vector[N] delta_gm;
    vector[N] delta_y;
    vector[N] delta_p; 
    vector[N_trt] log_beta_ref;
    
    // non-centered parameters
    log_p0 = log_p0_raw * sigma_log_p;
    log_p = log_p_raw * sigma_log_p;
    log_u = log_u_raw .* sigma_log_u[grp_p];
    
    // Exponents not vectorised in Stan
    delta_gm = exp(log(delta) * g_m);
    
    // Set reference class to zero
    log_beta_ref = append_row(rep_vector(0.0, N_grp), log_beta);
    
    // Use intercept for first obs.
    delta_y[t_1] = (1 - delta_gm[t_1]) .* 
                    (exp(log_p0[pop[t_1]] + log_u[plt[t_1]]) + obs_sd[obs[t_1]]);
    
    // Lagged obs. for remainder.
    delta_y[t] = (1 - delta_gm[t]) .* exp(log_y_std[t_m1]);
    
    // Remember predictions also offset
    delta_p = delta_gm .* (exp(log_p[pop] + 
                    log_beta_ref[trt] + log_u[plt]) + obs_sd[obs]);
    
    log_y_std ~ normal(log(delta_y + delta_p), 
                    (1 - delta_gm) / (1 - delta) .* sigma_log_e[grp]);
  }
  
  // priors
  log_p0_raw ~ normal(0, 1);
  log_p_raw ~ normal(0, 1);
  sigma_log_p ~ normal(0, 1);
  
  log_beta ~ normal(mu_log_b[grp_t], sigma_log_b[grp_t]);
  mu_log_b ~ normal(0, 1);
  sigma_log_b ~ normal(0, 1);
  
  delta ~ beta(4, 4);
  
  log_u_raw ~ normal(0, 1);
  sigma_log_u ~ normal(0, 1);
  sigma_log_e ~ normal(0, 1);
}
generated quantities{
  vector[T] y_hat[N_pop];
  
  for(i in 1:N_pop){
    y_hat[i, 1] = (1 - delta) * exp(log_p0_raw[i] * sigma_log_p);
    for(j in 2:T){
      y_hat[i, j] = (1 - delta) * y_hat[i, j - 1] + 
          delta * exp(log_p_raw[i] * sigma_log_p);
    }
  }
}
