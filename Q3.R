#Q3

#b)-----------------------------------------------------
log_post_boys <- function(gamma0, gamma1, gamma2, gamma3, tau, data) {
  if (tau <= 0) return(-Inf)
  
  t_obs <- c(0, 7, 14, 30*(1:6), 30*seq(12, 30, by=6))
  sigma <- sqrt(1 / tau)
  
  mu <- exp(gamma0) + 
        exp(-gamma1) * t_obs + 
        exp(gamma2) * (1 - exp(-exp(-gamma3 - gamma2) * t_obs))
  
  if (any(!is.finite(mu))) {
    return(-Inf)
  }
  
  log_lik <- 0
  for (j in 1:length(t_obs)) {
    log_lik <- log_lik + sum(dnorm(data[, j], mean = mu[j], sd = sigma, log = TRUE))
  }
  
  log_prior <- -log(tau)
  
  return(log_lik + log_prior)
}

log_post_boys(2,3,4,3,4,data)



#c)-----------------------------------------------------------------



neg_log_post_boys <- function(theta, data) {
  res <- log_post_boys(gamma0 = theta[1], 
                       gamma1 = theta[2], 
                       gamma2 = theta[3], 
                       gamma3 = theta[4], 
                       tau    = theta[5], 
                       data   = data)
  
  if (!is.finite(res)) return(1e9) 
  
  return(-res)
}

init_gamma0 <- log(beta0_boys)
init_gamma1 <- -log(beta1_boys)
init_gamma2 <- log(beta2_boys)
init_gamma3 <- -log(beta3_boys)
init_tau    <- 1 / (sigma_rough_boys^2)

theta_init <- c(init_gamma0, init_gamma1, init_gamma2, init_gamma3, init_tau)

fit_laplace <- nlm(f = neg_log_post_boys, 
                   p = theta_init, 
                   data = boys, 
                   hessian = TRUE)

theta_hat   <- fit_laplace$estimate
Sigma_theta <- solve(fit_laplace$hessian)

cat("--- Laplace Approximation Results ---\n")
cat("Theta_hat (Mode):\n")
print(round(theta_hat, 4))
cat("\nSigma_theta (Covariance Matrix):\n")
print(round(Sigma_theta, 6))