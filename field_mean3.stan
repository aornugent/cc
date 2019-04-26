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
  // Standardise
  vector[N] y_std = (y_obs ./ obs_sd[obs]) * obs_sd[12];
}
parameters{
  // intercepts
  vector<lower=0>[N_pop] p0_raw;
  
  // endpoints
  vector<lower=0>[N_pop] p_raw;
  real<lower=0> sigma_p;
  
  // predictors ordered so R- x N_grp at start
  vector<lower=0>[N_trt - N_grp] beta;
  vector<lower=0>[N_grp] mu_b;
  vector<lower=0>[N_grp] sigma_b;
  
  // correlation
  real<lower=0, upper=1> delta;
  
  // random effects
  vector<lower=0>[N_plt] u_raw;
  vector<lower=0>[N_grp] sigma_u;
  
  // obs. error
  vector<lower=0>[N_grp] sigma_e;
}
model{
  {
    vector[N_trt] beta_ref;
    vector[N] delta_y;
    vector[N] delta_p; 
    
    vector[N] delta_gm;
    vector[N_pop] p0;
    vector[N_pop] p;
    vector[N_plt] u;
    
    // Exponents not vectorised in Stan
    delta_gm = exp(log(delta) * g_m);
    
    // non-centered parameters
    p0 = p0_raw * sigma_p;
    p = p_raw * sigma_p;
    u = u_raw .* sigma_u[grp_p];
    
    // Set reference class to one
    beta_ref = append_row(rep_vector(1.0, N_grp), beta);
    
    // Use intercept for first obs.
    delta_y[t_1] = (1 - delta_gm[t_1]) .* (p0[pop[t_1]] .* u[plt[t_1]]);
    
    // Lagged obs. for remainder.
    delta_y[t] = (1 - delta_gm[t]) .* y_std[t_m1];
    
    // Multiplicative fixed and random effects
    delta_p = delta_gm .* (p[pop] .* u[plt] .* beta_ref[trt]);
    
    y_std ~ normal(delta_y + delta_p, 
                     (1 - delta_gm) / (1 - delta) .* sigma_e[grp]);
  }
  
  // priors
  p0_raw ~ normal(0, 1);
  p_raw ~ normal(0, 1);
  sigma_p ~ normal(0, 1);
  
  beta ~ normal(mu_b[grp_t], sigma_b[grp_t]);
  mu_b ~ normal(0, 1);
  sigma_b ~ normal(0, 1);
  
  delta ~ beta(4, 4);
  
  u_raw ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  sigma_e ~ normal(0, 1);
}
generated quantities{
  vector[T] y_hat[N_pop];
  
  for(i in 1:N_pop){
    y_hat[i, 1] = (1 - delta) * p0_raw[i] * sigma_p;
    for(j in 2:T){
      y_hat[i, j] = (1 - delta) * y_hat[i, j - 1] + delta * p_raw[i] * sigma_p;
    }
  }
}

