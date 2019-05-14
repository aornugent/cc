functions {
  // Eq. 14; Tjorve 2017. 
  vector gompertz_curve(vector x0, vector xK, 
                        vector rK, vector t, int n) {
                          
    // Power operation not vectorised in Stan
    vector[n] xt = log(xK) + log(x0 ./ xK) .* exp(-rK .* t);
    
    return(exp(xt));
  }
  
  // Eq. 7.38; Panik 2014.
  vector beta_curve(vector x0, vector xK, 
                    vector tK, vector tmax, 
                    vector t, int n) {
    
    vector[n] log_change = log(1 + (tK - t) ./ (tK - tmax)) + 
      log(t ./ tK) .* (tK ./ (tK - tmax));
      
    // back-transform and cap at 1
    vector[n] change;
    for(i in 1:n)
      change[i] = t[i] > tK[i] ? 1 : exp(log_change[i]);
    
    return(x0 + (xK - x0) .* change);
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
  // carrying capacities
  vector<lower=0>[N_pop] pK_raw;
  real<lower=0> sigma_p;
  
  // predictors ordered so R- x N_grp at start
  vector<lower=0>[N_trt - N_grp] beta_raw;
  vector<lower=0>[N_grp] mu_b;
  vector<lower=0>[N_grp] sigma_b;
  
  // correlation of error
  vector<lower=0, upper=1>[N_grp] delta;
  
  // plot random effects
  vector<lower=0>[N_plt] u_raw;
  vector<lower=0>[N_grp] sigma_u;
  
  // obs. error
  vector<lower=0>[N_grp] sigma_log_e;
}
transformed parameters {
    // non-centered parameters
    vector[N_pop] pK = pK_raw * sigma_p;
    vector[N_plt] u = u_raw .* sigma_u[grp_plt];
    vector[N_trt - N_grp] beta = mu_b[grp_trt] + beta_raw .* sigma_b[grp_trt];
}
model{
  {
    vector[N_trt] beta_ref;
    vector[N] mu_plt;
    vector[N] delta_gm;
    vector[N] log_mu;
    vector[N] delta_log_e;
    
    // Intercepts only
    // Set reference class to one
    beta_ref = append_row(rep_vector(1.0, N_grp), beta);
    
    // Multiplicative fixed and random effects
    mu_plt = pK[pop] .* u[plt] .* beta_ref[trt];
    

    // Stochastic model
    // Adjustment for irregular meas.
    delta_gm = exp(log(delta[grp]) .* gm);
    
    // Use mean for first meas.
    log_mu[m1] = log(delta_gm[m1] .* (pK[pop[m1]] .* u[plt[m1]]) + 
                    (1 - delta_gm[m1]) .* mu_plt[m1]);
    
    // Lagged meas. for reamainder
    log_mu[m] = log(delta_gm[m] .* y[m_m1] + (1 - delta_gm[m]) .* mu_plt[m]);


    // Measurement model
    // Uncertainty increases with time between meas.
    delta_log_e = (1 - delta_gm) ./ (1 - delta[grp]) .* sigma_log_e[grp];
    
    // Abundance follows lognormal distn.
    y[m_obs] ~ lognormal(log_mu[m_obs], delta_log_e[m_obs]);
    
    if(cen == 1){
      // Integrate out meas. below detection
      target += lognormal_lcdf(L[meas[m_mis]] | log_mu[m_mis], delta_log_e[m_mis]);
    }
  }
  
  // priors
  pK_raw ~ normal(0, 1);
  sigma_p ~ normal(0, 1);
  
  beta_raw ~ normal(0, 1);
  mu_b ~ normal(1, 1);
  sigma_b ~ normal(0, 1);
  
  delta ~ beta(1, 1);
  
  u_raw ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  sigma_log_e ~ normal(0, 1);
}