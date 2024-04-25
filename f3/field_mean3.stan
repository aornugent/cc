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
  int<lower=1> N_meas;  // measurment types
  int<lower=0> N_mis;   // censored
  
  // data
  vector<lower=1>[N] t;      // time
  vector<lower=0>[N] y;      // abundance
  vector<lower=0>[N] gm;     // gap between meas.
  vector<lower=0>[N_meas] L; // meas. limits
  int<lower=0, upper=1> cen; // switch for censoring
  int<lower=0, upper=1> model_variant; // varying capacity, rate or both

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
  vector<lower=0>[N_pop] p0_raw;
  
  // carrying capacity
  vector<lower=0>[N_pop] pK_raw;
  real<lower=0> sigma_p;

  // time that carrying capacity is reached
  vector<lower=0>[N_pop] log_tK_raw;
  real<lower=0> sigma_log_t;

  // time of max growth rate
  vector<lower=0, upper=1>[N_pop] tmax_raw;

  // predictors ordered so R- x N_grp at start
  vector<lower=0>[N_trt - N_grp] beta_raw;
  vector<lower=0>[N_grp] sigma_b;
  
  // correlation of error
  vector<lower=0, upper=1>[N_grp] delta;
  
  // plot random effectsQ
  vector<lower=0>[N_plt] u_raw;
  vector<lower=0>[N_grp] sigma_u;
  
  // obs. error
  vector<lower=0>[N_grp] sigma_log_e;
}
transformed parameters {
    // non-centered parameters
    vector[N_pop] p0 = p0_raw * sigma_p;
    vector[N_pop] pK = pK_raw * sigma_p;
    vector[N_pop] tK = T + exp(log_tK_raw * sigma_log_t);   // bound above T
    vector[N_pop] tmax = tmax_raw .* tK;                // bound below tK
    vector[N_plt] u = u_raw .* sigma_u[grp_plt];
    vector[N_trt - N_grp] beta = beta_raw .* sigma_b[grp_trt];
}
model{
  {
    vector[N_trt] beta_ref;
    vector[N] mu;
    vector[N] delta_gm;
    vector[N] lambda;
    vector[N] sigma_log_m;
    
    // Deterministic beta curve
    // Set reference class to one
    beta_ref = append_row(rep_vector(1.0, N_grp), beta);
    
     // Multiplicative fixed and random effects
    if(model_variant == 1) {
      mu = beta_curve(p0[pop] .* u[plt], 
                      pK[pop] .* u[plt] .* beta_ref[trt], 
                      tK[grp], tmax[grp], t, N);
    }
    else if (model_variant == 2){
      mu = beta_curve(p0[pop] .* u[plt], 
                      pK[pop] .* u[plt], 
                      tK[grp] .* beta_ref[trt],
                      tmax[grp] .* beta_ref[trt], t, N);
    }
    else if (model_variant == 3) {
      mu = beta_curve(p0[pop] .* u[plt], 
                      pK[pop] .* u[plt] .* beta_ref[trt], 
                      tK[grp] .* beta_ref[trt],
                      tmax[grp] .* beta_ref[trt], t, N);
    }
    
    // Stochastic model
    // Adjustment for irregular meas.
    delta_gm = exp(log(delta[grp]) .* gm);
    
    // Use latent obs for first meas.
    lambda[m1] = delta_gm[m1] .* p0[pop[m1]] .* u[plt[m1]] + 
                        (1 - delta_gm[m1]) .* mu[m1];
    
    // Lagged meas. for remainder
    lambda[m] = delta_gm[m] .* y[m_m1] + (1 - delta_gm[m]) .* mu[m];


    // Measurement model
    // Uncertainty increases with time between meas.
    sigma_log_m = (1 - delta_gm) ./ (1 - delta[grp]) .* sigma_log_e[grp];
    
    // Abundance follows lognormal distn.
    y[m_obs] ~ lognormal(log(lambda[m_obs]), sigma_log_m[m_obs]);
    
    if(cen == 1){
      // Integrate out meas. below detection
      target += lognormal_lcdf(L[meas[m_mis]] | log(lambda[m_mis]), sigma_log_m[m_mis]);
    }
  }
  
  // priors
  p0_raw ~ normal(1, 1);
  pK_raw ~ normal(1, 1);
  sigma_p ~ normal(0, 1);
  
  log_tK_raw ~ normal(0, 1);
  tmax_raw ~ beta(4, 4);
  sigma_log_t ~ cauchy(0, 1);
  
  beta_raw ~ normal(1, 1);
  sigma_b ~ normal(0, 1);
  
  delta ~ beta(1, 1);
  
  u_raw ~ normal(1, 1);
  sigma_u ~ normal(0, 1);
  sigma_log_e ~ normal(0, 1);
}
