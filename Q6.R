#Q6 
library(R2jags)
library(coda)
library(R2WinBUGS)
#taking data
w <- read.table("C:/Users/tsopa/Downloads/Weight.txt",header = TRUE)

# Measurement times in days
t.obs <- c(0, 7, 14, 30*(1:6), 30*seq(12, 30, by=6))


#keep 3 to 15 in order to not take the id and the fulo 
y_weights <- as.matrix(w[, 3:ncol(w)])


jags_data <- list(
  N=nrow(w),
  T = ncol(y_weights),
  t.obs= t.obs,         
  Male = w$Boy,
  y = y_weights 
)



mymodel=function(){
  for(i in 1:N){
    beta0[i]<-alpha0[i]+dbeta0.G*Male[i]
    beta1[i]<-beta1.F+dbeta1.G*Male[i]
    beta2[i]<-beta2.F+dbeta2.G*Male[i]
    beta3[i]<-beta3.F+dbeta3.G*Male[i]
    for (j in 1:T){
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


jags_params <- c("beta0.F", "beta1.F", "beta2.F", "beta3.F", 
                 "dbeta0.G", "dbeta1.G", "dbeta2.G", "dbeta3.G", 
                 "sigma", "sigma0.F", "dev")



set.seed(123)
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
print(jags_fit)


traceplot(jags_fit)

pdf("jags_traceplots_Q6.pdf", width = 8, height = 6)
traceplot(jags_fit, ask = FALSE)
dev.off()
