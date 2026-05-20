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
