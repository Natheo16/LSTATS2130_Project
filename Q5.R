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