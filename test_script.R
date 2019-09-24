library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("cc/run_model.R")

# Set field 
fields = 1:100
# fields = c(44, 53, 92)

# fit <- run_models("f3b", fields, censored = F, variant = 3)$fit
fit <- run_models("f2b", fields, censored = T, variant = 1, chains = 4, iter = 600)$fit
