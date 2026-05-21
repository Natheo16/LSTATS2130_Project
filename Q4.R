#Q4

library(MASS)
library(coda)

set.seed(123)

M <- 10000
burnin <- 2000
N_obs <- nrow(boys) * ncol(boys) 

gamma_post <- matrix(NA, nrow = M, ncol = 4)
tau_post   <- numeric(M)

gamma_post[1, ] <- theta_hat[1:4]
tau_post[1]     <- theta_hat[5]

c_tune     <- (2.38^2) / 4 
Sigma_prop <- c_tune * Sigma_theta[1:4, 1:4]

n_accept <- 0

log_lik_gamma <- function(gamma, tau, data, t_obs) {
  if (any(!is.finite(gamma))) return(-Inf)
  
  sigma <- sqrt(1/tau)
  mu <- exp(gamma[1]) + exp(-gamma[2]) * t_obs + exp(gamma[3]) * (1 - exp(-exp(-gamma[4] - gamma[3]) * t_obs))
  
  if (any(!is.finite(mu))) return(-Inf)
  
  ll <- 0
  for (j in 1:length(t_obs)) {
    ll <- ll + sum(dnorm(data[, j], mean = mu[j], sd = sigma, log = TRUE))
  }
  
  return(ll)
}

#  MCMC Loop 
for (m in 2:M) {
  
  gamma_prev <- gamma_post[m-1, ]
  tau_prev   <- tau_post[m-1]
  
  gamma_prop <- mvrnorm(1, mu = gamma_prev, Sigma = Sigma_prop)
  
  log_r <- log_lik_gamma(gamma_prop, tau_prev, boys, t.obs) - 
           log_lik_gamma(gamma_prev, tau_prev, boys, t.obs)
  
  if (is.na(log_r)) log_r <- -Inf
  
  if (runif(1) <= exp(log_r)) {
    gamma_post[m, ] <- gamma_prop
    n_accept <- n_accept + 1
  } else {
    gamma_post[m, ] <- gamma_prev
  }
  
  gamma_curr <- gamma_post[m, ]
  mu_curr <- exp(gamma_curr[1]) + exp(-gamma_curr[2]) * t.obs + 
             exp(gamma_curr[3]) * (1 - exp(-exp(-gamma_curr[4] - gamma_curr[3]) * t.obs))
  
  SSR <- 0
  for (j in 1:length(t.obs)) {
    SSR <- SSR + sum((boys[, j] - mu_curr[j])^2)
  }
  
  tau_post[m] <- rgamma(1, shape = N_obs / 2, rate = SSR / 2)
}

#  ii and iii. Diagnostics and ESS 
gamma_samples <- gamma_post[-(1:burnin), ]
tau_samples   <- tau_post[-(1:burnin)]

chain <- mcmc(cbind(gamma_samples, tau_samples))
colnames(chain) <- c("gamma0", "gamma1", "gamma2", "gamma3", "tau")

cat("Acceptance rate for gamma:", n_accept / M, "\n\n")

plot(chain)
print(geweke.diag(chain))

cat("\nEffective Sample Size:\n")
print(round(effectiveSize(chain), 0))

#  iv. Inference and Back-transformation 
beta_chains <- cbind(
  beta0 = exp(gamma_samples[, 1]),
  beta1 = exp(-gamma_samples[, 2]),
  beta2 = exp(gamma_samples[, 3]),
  beta3 = exp(-gamma_samples[, 4])
)

results <- apply(beta_chains, 2, function(x) {
  c(Point_Est = median(x), 
    Lower_95  = quantile(x, 0.025), 
    Upper_95  = quantile(x, 0.975))
})

cat("\n--- Point Estimates and 95% CI for Beta parameters ---\n")
print(round(t(results), 4))
pdf("traceplots.pdf", width = 8, height = 6)
plot(chain)
dev.off()