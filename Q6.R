#Q6 
library(R2jags)
library(coda)
library(R2WinBUGS)
#taking data
w <- read.table("C:/Users/tsopa/Downloads/Weight.txt",header = TRUE)

# Measurement times in days
t.obs <- c(0, 7, 14, 30*(1:6), 30*seq(12, 30, by=6))


# Extract weight measurements (columns 3 to 15) into a matrix
y_weights <- as.matrix(w[, 3:ncol(w)])

# Prepare data list for JAGS
jags_data <- list(
  N=nrow(w),
  T = ncol(y_weights),
  t.obs= t.obs,         
  Male = w$Boy,
  y = y_weights 
)



mymodel=function(){
  for(i in 1:N){
    # Gender-specific parameters (Random effect alpha0 applied only to intercept)
    beta0[i]<-alpha0[i]+dbeta0.G*Male[i]
    beta1[i]<-beta1.F+dbeta1.G*Male[i]
    beta2[i]<-beta2.F+dbeta2.G*Male[i]
    beta3[i]<-beta3.F+dbeta3.G*Male[i]
    for (j in 1:T){
      # Non-linear exponential growth curve equation
      mu[i,j]<-beta0[i]+beta1[i]*t.obs[j]+beta2[i]*(1-exp(-beta3[i]/beta2[i]*t.obs[j]))
      y[i,j]~dnorm(mu[i,j],tau)
      #likehood contribution
      Lik[i,j] <- dnorm(y[i,j],mu[i,j],tau)
    }
    alpha0[i] ~ dnorm(beta0.F,tau0.F) ## Varying intercepts
  }
  #Priors
  beta0.F ~ dnorm(0,0.0001) ; tau0.F ~ dgamma(0.0001,0.0001)
  beta1.F ~ dnorm(0,0.0001)
  beta2.F ~ dnorm(0,0.0001)
  beta3.F ~ dnorm(0,0.0001)
  dbeta0.G ~ dnorm(0,0.0001) ; 
  dbeta1.G ~ dnorm(0,0.0001);
  dbeta2.G ~ dnorm(0,0.0001) ; 
  dbeta3.G ~ dnorm(0,0.0001)
  tau ~ dgamma(0.0001,0.0001)
  sigma0.F <- 1/sqrt(tau0.F)
  sigma <- 1/sqrt(tau)
  ## Deviance
  dev <- -2*sum(log(Lik))
}

# Write the model function into a text file for JAGS
model.file = "model_exercise6.bug"
write.model(mymodel, model.file)

#function putting the values we know from previous exercises
jags_inits <- function() {
  list(
    beta0.F = rnorm(1, mean = 2.95, sd = 0.1),
    beta1.F = rnorm(1, mean = 0.005, sd = 0.001),
    beta2.F = rnorm(1, mean = 5.3, sd = 0.1),
    beta3.F = rnorm(1, mean = 0.037, sd = 0.001),
    dbeta0.G = rnorm(1, mean = 0, sd = 0.1),
    dbeta1.G = rnorm(1, mean = 0, sd = 0.001),
    dbeta2.G = rnorm(1, mean = 0, sd = 0.1),
    dbeta3.G = rnorm(1, mean = 0, sd = 0.001),
    tau = runif(1, min = 5, max = 15),
    tau0.F = runif(1, min = 5, max = 15),
    alpha0 = rnorm(nrow(w), mean = 2.95, sd = 0.1) 
  )
}

# Target parameters to monitor and save from the posterior distribution
jags_params <- c("beta0.F", "beta1.F", "beta2.F", "beta3.F", 
                 "dbeta0.G", "dbeta1.G", "dbeta2.G", "dbeta3.G", 
                 "sigma", "sigma0.F", "dev")



set.seed(123)
# Run the JAGS simulation using 3 parallel Markov chains
jags_fit <- jags(
  data = jags_data,
  inits = jags_inits,
  parameters.to.save = jags_params,
  model.file = model.file,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 1
)  
# Display posterior summaries table (Means, SDs, CIs, Rhat, DIC)
print(jags_fit)


traceplot(jags_fit)
# Save traceplots into a PDF file for reporting
pdf("jags_traceplots_Q6.pdf", width = 8, height = 6)
traceplot(jags_fit, ask = FALSE)
dev.off()








# ==============================================================================
# QUESTION 6(d): Posterior Predictive Distribution for a Boy Aged 90 Days
# ==============================================================================

# 1. Define the updated model function including the 90-day predictive logic
mymodel_d = function(){
  for(i in 1:N){
    beta0[i] <- alpha0[i] + dbeta0.G * Male[i]
    beta1[i] <- beta1.F + dbeta1.G * Male[i]
    beta2[i] <- beta2.F + dbeta2.G * Male[i]
    beta3[i] <- beta3.F + dbeta3.G * Male[i]
    
    for (j in 1:T){
      mu[i,j] <- beta0[i] + beta1[i] * t.obs[j] + beta2[i] * (1 - exp(-beta3[i] / beta2[i] * t.obs[j]))
      y[i,j] ~ dnorm(mu[i,j], tau)
    }
    alpha0[i] ~ dnorm(beta0.F, tau0.F) ## Varying intercepts
  }
  
  # --- Posterior Predictive Logic for a NEW Boy Aged 90 Days ---
  # Draw a new random intercept for an unobserved child from the population
  alpha0.new.Boy ~ dnorm(beta0.F, tau0.F)
  
  # Compute the structural curves for a boy (Male = 1) at t = 90 days
  beta0.new.Boy <- alpha0.new.Boy + dbeta0.G * 1
  beta1.new.Boy <- beta1.F + dbeta1.G * 1
  beta2.new.Boy <- beta2.F + dbeta2.G * 1
  beta3.new.Boy <- beta3.F + dbeta3.G * 1
  
  # Expected mean weight at 90 days
  mu.boy.90 <- beta0.new.Boy + beta1.new.Boy * 90 + beta2.new.Boy * (1 - exp(-beta3.new.Boy / beta2.new.Boy * 90))
  
  # Final simulated weight incorporating measurement error
  weight.boy.90 ~ dnorm(mu.boy.90, tau)
  
  # --- Priors ---
  beta0.F ~ dnorm(0, 0.0001)
  tau0.F ~ dgamma(0.0001, 0.0001)
  beta1.F ~ dnorm(0, 0.0001)
  beta2.F ~ dnorm(0, 0.0001)
  beta3.F ~ dnorm(0, 0.0001)
  
  dbeta0.G ~ dnorm(0, 0.0001)
  dbeta1.G ~ dnorm(0, 0.0001)
  dbeta2.G ~ dnorm(0, 0.0001)
  dbeta3.G ~ dnorm(0, 0.0001)
  
  tau ~ dgamma(0.0001, 0.0001)
  sigma0.F <- 1 / sqrt(tau0.F)
  sigma <- 1 / sqrt(tau)
}

# 2. Write the updated function to the bug file
model.file.d = "model_exercise6_d.bug"
write.model(mymodel_d, model.file.d)

# 3. Update parameters to track the new predictive variable
jags_params_d <- c("beta0.F", "beta1.F", "beta2.F", "beta3.F", 
                   "dbeta0.G", "dbeta1.G", "dbeta2.G", "dbeta3.G", 
                   "sigma", "sigma0.F", "weight.boy.90")

# 4. Run the JAGS simulation using your verified jags_data and jags_inits
set.seed(123)
jags_fit_d <- jags(
  data = jags_data,
  inits = jags_inits,
  parameters.to.save = jags_params_d,
  model.file = model.file.d,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 1
)

# 5. Print the summary specifically for the 90-day predictive weight
print("--- POSTERIOR PREDICTIVE MEASURES FOR A BOY AT 90 DAYS ---")
print(jags_fit_d$BUGSoutput$summary["weight.boy.90", c("mean", "sd", "2.5%", "50%", "97.5%")])

# 6. Visualize the Predictive Distribution (Academic Black & White Density Plot)
mcmc_samples_90 <- jags_fit_d$BUGSoutput$sims.list$weight.boy.90
boy_stats <- jags_fit_d$BUGSoutput$summary["weight.boy.90", ]

plot(density(mcmc_samples_90), 
     main = "Posterior Predictive Distribution for a New Boy at 90 Days",
     xlab = "Predicted Weight (kg)", ylab = "Density", 
     col = "black", lwd = 2) 

# Add line markers for Median and 95% Credible Intervals using styles instead of colors
abline(v = boy_stats["50%"], col = "black", lty = 1, lwd = 2) 
abline(v = boy_stats["2.5%"], col = "black", lty = 2, lwd = 1.5) 
abline(v = boy_stats["97.5%"], col = "black", lty = 2, lwd = 1.5) 

legend("topright", legend = c("Predictive Density", "Median", "95% CI Bounds"),
       col = c("black", "black", "black"), lty = c(1, 1, 2), lwd = c(2, 2, 1.5))
