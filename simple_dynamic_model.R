
library(tidyverse)
library(rstan)

# load Gompertz function
expose_stan_functions("cc/f2/field_mean2.stan")

# Simple dynamic model testing auto-correlatoin
n = 100
x0 = 6
xK = 10
delta = .9

x = accumulate(rep(delta, n), 
       ~ rlnorm(1, log(.y * .x + (1 - .y) * xK), log(1.01)), 
       .init = x0)

plot(0:n, x, col = "blue")

# Irregular obs
m = sort(sample(1:n, 40))
g_m = m - lag(m, default = 0)

# Adjust variance for time between obs, can sometimes be much larger
x_m = accumulate(delta^g_m, 
       ~ rlnorm(1, log(.y * .x + (1 - .y) * xK), 
              ((1 - delta^g_m) / (1 - delta)) * 0.01), 
       .init = x0)

points(c(0, m), x_m, col = "red")

# Deterministic Gompertz function
# Eqn 14: Tjorve (2017)
rK = .2

x_g = accumulate(1:n, ~ xK * (x0 / xK)^exp(-rK * .y), .init = x0)
x_g = accumulate(1:n, ~ gompertz_curve(x0, xK, rK, .y, 1), .init = x0)

points(0:n, x_g, col = "black")