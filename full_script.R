# Q1

library(ggplot2)
library(tidyr)
library(dplyr)

data <- read.table("Weight.txt", header = TRUE)

# Measurement times in days
t.obs <- c(0, 7, 14, 30*(1:6), 30*seq(12, 30, by=6))

# Rename columns to actual time values
colnames(data) <- c("Boy", t.obs)

# Add subject ID
data$id <- 1:nrow(data)

# Reshape to usable format
data_long <- data %>%
  pivot_longer(cols = -c(id, Boy), names_to = "time", values_to = "weight") %>%
  mutate(time = as.numeric(time),
         sex  = ifelse(Boy == 1, "Boys", "Girls"))

# --- Plot 1: Individual trajectories ---
p1 <- ggplot(data_long, aes(x = time, y = weight, group = id)) +
  geom_line(alpha = 0.3, linewidth = 0.4) +
  facet_wrap(~sex) +
  labs(x = "Age (days)", y = "Weight (kg)",
       title = "Individual weight trajectories by sex") +
  theme_bw()

ggsave("Q1_trajectories.png", p1, width = 10, height = 5, dpi = 150)

# --- Plot 2: Mean weight over time ---
mean_data <- data_long %>%
  group_by(sex, time) %>%
  summarise(mean_weight = mean(weight), .groups = "drop")

p2 <- ggplot(mean_data, aes(x = time, y = mean_weight, color = sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Boys" = "steelblue", "Girls" = "tomato")) +
  labs(x = "Age (days)", y = "Mean weight (kg)",
       title = "Mean weight over time by sex", color = "Sex") +
  theme_bw()

ggsave("Q1_mean_growth.png", p2, width = 7, height = 5, dpi = 150)



# Q2

# --- 2d: Rough estimates from the graphs ---

data <- read.table("Weight.txt", header = TRUE)
t.obs <- c(0, 7, 14, 30*(1:6), 30*seq(12, 30, by=6))

boys  <- data[data$Boy == 1, -1]
girls <- data[data$Boy == 0, -1]

mean_boys  <- colMeans(boys)
mean_girls <- colMeans(girls)


beta0_boys <- mean_boys[1]
beta1_boys <- (mean_boys[13] - mean_boys[11]) / (t.obs[13] - t.obs[11])
init_slope_boys <- (mean_boys[2] - mean_boys[1]) / (t.obs[2] - t.obs[1])
beta3_boys <- init_slope_boys - beta1_boys
beta2_boys <- mean_boys[11] - beta0_boys - beta1_boys * t.obs[11]

cat("--- Boys rough estimates ---\n")
cat(sprintf("beta0 ≈ %.3f kg\n",    beta0_boys))
cat(sprintf("beta1 ≈ %.5f kg/day\n", beta1_boys))
cat(sprintf("beta2 ≈ %.3f kg\n",    beta2_boys))
cat(sprintf("beta3 ≈ %.5f kg/day\n", beta3_boys))

# Girls:
beta0_girls <- mean_girls[1]
beta1_girls <- (mean_girls[13] - mean_girls[11]) / (t.obs[13] - t.obs[11])
init_slope_girls <- (mean_girls[2] - mean_girls[1]) / (t.obs[2] - t.obs[1])
beta3_girls <- init_slope_girls - beta1_girls
beta2_girls <- mean_girls[11] - beta0_girls - beta1_girls * t.obs[11]

cat("\n--- Girls rough estimates ---\n")
cat(sprintf("beta0 ≈ %.3f kg\n",    beta0_girls))
cat(sprintf("beta1 ≈ %.5f kg/day\n", beta1_girls))
cat(sprintf("beta2 ≈ %.3f kg\n",    beta2_girls))
cat(sprintf("beta3 ≈ %.5f kg/day\n", beta3_girls))

# --- 2e: Rough estimate of sigma ---

sigma_rough_boys  <- sd(boys[, 10])   # 12m column
sigma_rough_girls <- sd(girls[, 10])

cat(sprintf("\nSample SD at 12 months (boys):  %.3f kg\n", sigma_rough_boys))
cat(sprintf("Sample SD at 12 months (girls): %.3f kg\n", sigma_rough_girls))





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

#Q5
library(R2jags)

jags_data <- list(
  boys = boys,           
  t_obs = t.obs,         
  N_boys = nrow(boys),
  N_time = length(t.obs)
)

jags_inits <- function() {
  list(
    beta0 = rnorm(1, mean = 2.95, sd = 0.1),
    beta1 = rnorm(1, mean = 0.005, sd = 0.001),
    beta2 = rnorm(1, mean = 5.3, sd = 0.1),
    beta3 = rnorm(1, mean = 0.037, sd = 0.001),
    tau   = runif(1, min = 5, max = 15)
  )
}

jags_params <- c("beta0", "beta1", "beta2", "beta3", "sigma")

set.seed(123)
jags_fit <- jags(
  data = jags_data,
  inits = jags_inits,
  parameters.to.save = jags_params,
  model.file = "model_boys.bug",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 1
)

print(jags_fit)

traceplot(jags_fit)

pdf("jags_traceplots.pdf", width = 8, height = 6)
traceplot(jags_fit, ask = FALSE)
dev.off()

