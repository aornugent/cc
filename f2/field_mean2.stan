functions {
  // Eq. 14; Tjorve 2017. 
  vector gompertz_curve(vector x0, vector xK, 
                        vector rK, vector t, int n) {
                          
    // Power operation not vectorised in Stan
    vector[n] xt = log(xK) + log(x0 ./ xK) .* exp(-rK .* t);
    
    return(exp(xt));
  }
}
data {
  // sizes
  int<lower=1> T;       // timespan
  int<lower=1> N;       // total
  int<lower=1> N_grp;   // groups
  int<lower=1> N_trt;   // contrasts
  int<lower=1> N_plt;   // experiments
  int<lower=1> N_pop;   // populations
  int<lower=1> N_meas;  // meas. types
  int<lower=0> N_mis;   // censored
  
  // data
  vector<lower=1>[N] t;      // time
  vector<lower=0>[N] y;      // abundance
  vector<lower=0>[N] gm;     // gap between meas.
  vector<lower=0>[N_meas] L; // limits
  int<lower=0, upper=1> cen; // switch for censoring

  // meas. indices
  int<lower=1, upper=N> m1[N_plt];
  int<lower=1, upper=N> m[N - N_plt];
  int<lower=1, upper=N> m_m1[N - N_plt];
  int<lower=1, upper=N> m_obs[N - N_mis];
  int<lower=0, upper=N> m_mis[N_mis];
  
  // meas., group, treatment, plot & pop. indices
  int<lower=1, upper=N_meas> meas[N];
  int<lower=1, upper=N_grp> grp[N];
  int<lower=1, upper=N_trt> trt[N];
  int<lower=1, upper=N_grp> grp_trt[N_trt - N_grp];
  int<lower=1, upper=N_plt> plt[N];
  int<lower=1, upper=N_grp> grp_plt[N_plt];
  int<lower=1, upper=N_pop> pop[N]; 
  int<lower=1, upper=N_grp> grp_pop[N_pop];
}
parameters{
  // latent starting point
  vector<lower=0>[N_pop] p0;
  vector<lower=0>[N_grp] mu_p0;
  
  // mean abundance
  vector<lower=0>[N_pop] pK;
  vector<lower=0>[N_grp] mu_pK;
  real<lower=0> sigma_p;
  
  // vector<lower=0>[N_grp] rK;
  vector<lower=0>[N_pop] rK;
  vector<lower=0>[N_grp] mu_rK;
  real<lower=0> sigma_r;
  
  // predictors ordered so R- x N_grp at start
  vector<lower=0>[N_trt - N_grp] beta;
  real<lower=0> sigma_b;
  
  // correlation of error
  vector<lower=0, upper=1>[N_grp] delta;
  
  // plot random effectsQ
  vector<lower=0>[N_plt] u;
  real<lower=0> sigma_u;
  
  // obs. error
  vector<lower=0>[N_pop] sigma_e;
  vector<lower=0>[N_grp] mu_sigma_e;
  real<lower=0> sigma_sigma_e;
}
model{
  {
    vector[N_trt] beta_ref;
    vector[N] lambda;
    vector[N] mu;
    vector[N] delta_gm;
    vector[N] log_mu;
    vector[N] sigma_m;
    vector[N] log_sigma_m_sq;
    
    // Intercepts only
    // Set reference class to one
    beta_ref = append_row(rep_vector(1.0, N_grp), beta);
    
    // Multiplicative fixed and random effects
    lambda = gompertz_curve(p0[pop] .* u[plt], 
                            pK[pop] .* u[plt] .* beta_ref[trt], 
                            rK[pop], t, N);

    // Stochastic model
    // Adjustment for irregular meas.
    delta_gm = exp(log(delta[grp]) .* gm);
    
    // Use latent obs for first meas.
    mu[m1] = lambda[m1] + delta_gm[m1] .*(p0[pop[m1]] .* u[plt[m1]] - lambda[m1]);
    
    // Lagged meas. for remainder
    mu[m] = lambda[m] + delta_gm[m] .* (y[m_m1] - lambda[m]);
    

    // Measurement model
    // Uncertainty increases with time between meas.
    log_sigma_m_sq = log(1 + ((1 - delta_gm) ./ (1 - delta[grp]) .* 
                                      square(sigma_e[pop])) ./ square(mu));
    
   // Transform to log scale
    log_mu = log(mu) - 0.5 * log_sigma_m_sq;
    
    // Abundance follows lognormal distn.
    y[m_obs] ~ lognormal(log_mu[m_obs], sqrt(log_sigma_m_sq[m_obs]));
    
    if(cen == 1){
      // Integrate out meas. below detection
      target += lognormal_lcdf(L[meas[m_mis]] | 
                                log_mu[m_mis],
                                sqrt(log_sigma_m_sq[m_mis]));
    }
  }
  
  // priors
  p0 ~ normal(mu_p0[grp_pop], sigma_p);
  pK ~ normal(mu_pK[grp_pop], sigma_p);
  mu_p0 ~ normal(1, 1);
  mu_pK ~ normal(1, 1);
  sigma_p ~ normal(0, 1);
  
  // rK ~ normal(0, 1);
  rK ~ normal(mu_rK[grp_pop], sigma_r);
  mu_rK ~ normal(0, 1);
  sigma_r ~ normal(0, 1);

  beta ~ normal(1, sigma_b);
  sigma_b ~ normal(0, 1);
  
  delta ~ beta(1, 1);
  
  u ~ normal(1, sigma_u);
  sigma_u ~ normal(0, 1);
  sigma_e ~ normal(mu_sigma_e[grp_pop], sigma_sigma_e);
  mu_sigma_e ~ normal(0, 1);
  sigma_sigma_e ~ normal(0, 1);
}
