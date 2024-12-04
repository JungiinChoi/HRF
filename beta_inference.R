# Inference for beta for different cases

# Load Rscripts
setwd("Pain")
source("../Fit_Canonical_HRF.R")
source("../Fit_Spline.R")
source("../get_parameters2.R")
source("../spm_hrf.R")
source("../Det_Logit.R")
source("../tor_make_deconv_mtx3.R")

# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(sp)
library(gstat)
library(raster)
library(stats)
library(mgcv)  # for thin plate regression spline
library(akima)
library(Matrix)

# Simulation: DD + one beta


# Model: Canonical DD (two betas)

# I. 2D Square Region

# Define the test function f(x, z)
f_function_1 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (3.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.2)^2 / sigma_x^2) - ((z - 0.3)^2 / sigma_z^2))
  term2 <- (2.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.7)^2 / sigma_x^2) - ((z - 0.8)^2 / sigma_z^2))
  return(term1 + term2)
}

# Set parameters
sigma_x <- 0.3
sigma_z <- 0.4
n_points <- 500  # number of points to generate

# Generate random (x, z) points within the unit square
set.seed(123)
grid <- expand.grid(x = seq(0, 1, length.out = sqrt(n_points)), y = seq(0, 1, length.out = sqrt(n_points)))

x <- grid[,1]
z <- grid[,2]

# Create a data frame from the vectors
data <- data.frame(x1 = x, x2 = z)

# Use ggplot to plot the points
ggplot(data, aes(x = x1, y = x2)) +
  geom_point(color = "blue") +  # Plot points
  labs(title = "Voxels", x = "X1", y = "X2") +  # Add labels
  theme_minimal()  # Use a clean 

# Create a grid for plotting the true function
grid_size <- 50
x_grid <- seq(0, 1, length.out = grid_size)
z_grid <- seq(0, 1, length.out = grid_size)
beta1_grid <- outer(x_grid, z_grid, function(x, z) f_function_1(x, z, sigma_x, sigma_z))

# Draw contour plots'
# Contour plot of the true function
contour(x_grid, z_grid, beta1_grid, main = "True Beta1", xlab = "x1", ylab = "x2")

# Case I: Large noise (fmri data has large noise!)
# I. independent(non-spatial) model

# Settings
TR = 1
T = 40
t = seq(1, T, by = TR)    # samples at which to get Logit HRF Estimate
FWHM = 4                  # FWHM for residual scan
pval = 0.01
df = 600
alpha = 0.001
len = 640

# function generating error 
noise_arp <- function(n, phi = c(0.5, 0.1), sigma = 1, e = rnorm(n, 0, sigma)){
  p <- length(phi)
  w <- numeric(n)
  
  for (i in 1:n) {
    if (i <= p) {
      w[i] <- e[i]
    } else {
      arpart <- sum(phi * w[(i-1):(i-p)], na.rm = TRUE)
      w[i] <- arpart + e[i]
    }
  }
  return(w)
}

# Plot the hrf
beta_true_1 <- seq(0.317,9.044,length = 30)
hrfdf <- data.frame()
hrfmat <- matrix(nrow = length(beta_true_1), ncol = 40)
for (i in 1:length(beta_true_1)){
  true_hrf <- CanonicalBasisSet(TR)$h * beta_true_1[i]
  hrfdf <- rbind(hrfdf, data.frame(index = rep(i,40), x = 0:39, y = true_hrf))
}


hrfdf$index <- factor(hrfdf$index)

ggplot(hrfdf, aes(x, y, group = index, colour = index)) +
  labs(title = "HRF for Different Betas", x = "t", y = "y") +
  geom_line() 

# Onset
b <- 1
R <- c(13, 14, 29, 44, 107, 125, 160, 171, 174, 190, 191, 
       206, 215, 232, 237, 262, 277, 292, 296, 346, 354, 
       367, 375, 382, 384, 398, 409, 462, 469, 475, 501, 
       520, 527, 566, 577, 629)
Run <- rep(0, 640)
Run[R] <- 1
Runc <- list(Run)
onset <- rep(-100,640)
onset[R] <- R

# Simulations
set.seed(906)
tc_mg <- c()

sim_k <- 10
error_sig <- c(5,10,30,50,100)
beta_DD <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1))
sig_hat <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1))

for (j in 1:length(error_sig)){
  for (k in 1:length(beta1)){
    true_sig <-  beta1[k] * conv(Run, CanonicalBasisSet(TR)$h)[1:len]
    xsecs <- 0:40
    tc_mat <- matrix(0,nrow = len, ncol = sim_k)
    
    for (i in 1:sim_k){
      tc_noise <- noise_arp(n = len, phi = c(0.3, 0), sigma = error_sig[j])
      tc <- true_sig + tc_noise
      if (i ==1){
        tc_mg <- c(tc_mg, tc)
      }
      tc_mat[,i] <- tc
      
      #Canonical HRF
      Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Runc, T, 1)
      sig_hat[k,i,j,] <- sqrt(sum(Canonical_HRF_fitted$e^2)/(length(Canonical_HRF_fitted$e)-1))
      beta_DD[k,i,j,] <- Canonical_HRF_fitted$info$b[1,]
    }
    print(k)
  }
}


ratio <- rep(0,length(beta1))

for (i in 1:length(beta1)){
  mean_value <- mean(beta_DD[i,,1,1])
  std_error <- sd(beta_DD[i,,1,1])
  t_value <- qnorm(0.975)
  
  # Calculate the 95% confidence interval
  ci_lower <- mean_value - t_value * std_error
  ci_upper <- mean_value + t_value * std_error
  
  if (beta1[i] >= ci_lower & beta1[i] <= ci_upper){
    ratio[i] <- 1
  }
}

mean(ratio)

# Inference on sigma_k

# Calculate the sample mean
mean_value <- mean(sig_hat)

# Calculate the standard error of the mean
n <- length(sig_hat)
std_error <- sd(sig_hat)

# Get the critical t-value for a 95% confidence level
# Using n-1 degrees of freedom
t_value <- qt(0.975, df = n - 1)

# Calculate the 95% confidence interval
ci_lower <- mean_value - t_value * std_error
ci_upper <- mean_value + t_value * std_error

# Print the results
mean_value
cat("95% Confidence Interval: [", ci_lower, ", ", ci_upper, "]\n")



# Visualize timecourse

# Calculate mean and 95% confidence interval for each time point
time_points <- 1:nrow(tc_mat)
means <- apply(tc_mat, 1, mean)
sds <- apply(tc_mat, 1, sd)
n <- ncol(tc_mat)

# Calculate standard error and confidence intervals
stderr <- sds / sqrt(n)
conf_interval <- qt(0.975, df = n - 1) * stderr

# Create data frame for ggplot
data <- data.frame(
  Time = time_points,
  Mean = means,
  CI_Lower = means - conf_interval,
  CI_Upper = means + conf_interval
)

# Plot using ggplot
ggplot(data, aes(x = Time, y = Mean)) +
  geom_line(color = "blue", linewidth = 0.7) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, fill = "blue") +
  labs(title = expression("Timecourse Data for " ~ sigma == 50),
       x = "t",
       y = "y") +
  theme_minimal()



# Estimate beta independently

ind_beta1 <- apply(beta_DD[,,1],1,mean)
ind_dist <- (beta1 - ind_beta1)

df_fitted <- data.frame(x1 = x, x2 = z, ind_beta1 = ind_beta1)

#create plot

fld <- with(df_ind, interp(x = x1, y = x2, z = ind_beta1))
fld$x <- c(fld$x, 1.000001)
fld$y <- c(fld$y, 1.000001)
fld$z <- rbind(cbind(fld$z, rep(10, 40)), rep(10,41))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "Independent Beta",
               key.title = title(main = "beta1 ", cex.main = 1))

contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(ind_beta1, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
        main = "Independent Beta", xlab = "x1", ylab = "x2")



# II. Gaussian Kernel Smoothing + using smoothing var

# Define bandwidth for the Gaussian kernel
bandwidth <- 0.045

# Sample data: Replace with your actual beta values and coordinates
x_coords <- x
y_coords <- z
beta_values <- ind_beta1

# Combine into a data frame
data <- data.frame(x = x_coords, y = y_coords)

data <- data.frame(x = x_coords, y = y_coords, beta = beta_values)

# Initialize a matrix to store smoothed beta values
smoothed_beta <- rep(0,length(x))

# Function to compute Gaussian kernel weight
gaussian_weight <- function(dx, dy, bandwidth) {
  exp(- (dx^2 + dy^2) / (2 * bandwidth^2))
}


# save timecourse data
set.seed(906)
tc_mg <- c()

sim_k <- 100
error_sig <- c(5,10,30,50,100)

tc_mat <- array(0, dim = c(len, length(beta1), sim_k, length(error_sig)))

for (j in 1:length(error_sig)){
  for (k in 1:length(beta1)){
    true_sig <-  beta1[k] * conv(Run, CanonicalBasisSet(TR)$h)[1:len]
    xsecs <- 0:40
    for (i in 1:sim_k){
      tc_noise <- noise_arp(n = len, phi = c(0.3, 0), sigma = error_sig[j])
      tc <- true_sig + tc_noise
      tc_mat[,k,i,j] <- tc
    }
    print(k)
  }
}


# Calculate smoothed beta and Y values

smoothed_beta <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1))
smoothed_Y <- array(0, dim = c(len, length(beta1), sim_k, length(error_sig)))
weights_total <- matrix(0, length(beta1), length(beta1))

for (j in 1:length(error_sig)){
  for (k in 1:length(beta1)){
    # Get the current grid point
    gx <- x[k]
    gy <- z[k]
    # Calculate weights based on the Gaussian kernel
    weights <- mapply(function(x, y) gaussian_weight(gx - x, gy - y, bandwidth), data$x, data$y)
    weights_total[k,] <- weights
    
    for (i in 1:sim_k){
      # Compute weighted average of beta values
      smoothed_beta[k,i,j,1] <- sum(weights * beta_DD[,i,j,1]) / sum(weights)
      for (l in 1:len){
        smoothed_Y[l,k,i,j] <- sum(weights * tc_mat[l,,i,j]) / sum(weights)
      }
    }
    print(k)
  }
}

var_weight <- apply(weights_total, 1, function(x) sqrt(sum(x^2)/(sum(x))^2))
contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(var_weight, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
        main = expression(sqrt(frac(sum(K[i]^2), (sum(K[i]))^2))), xlab = "x1", ylab = "x2", nlevels = 10)
hist(var_weight, main="")
title(main = expression(sqrt(frac(sum(K[i]^2), (sum(K[i]))^2))))


sig_hat_tilde <- array(0, dim = c(length(error_sig), length(beta1), sim_k))
X <- Canonical_HRF_fitted$info$X[,2]

for (i in 1:length(error_sig)){
  for (j in 1:length(beta1)){
    for (k in 1:sim_k){
      sig_hat_tilde[i,j,k] <- sum((smoothed_Y[,j,k,i] - X * smoothed_beta[j,k,i,1])^2)/(len - 1)
    }
  }
  
  mean_value <- mean(sig_hat_tilde[i,,])
  std_error <- sd(sig_hat_tilde[i,,])
  t_value <- qnorm(0.975)
  
  # Calculate the 95% confidence interval
  ci_lower <- mean_value - t_value * std_error
  ci_upper <- mean_value + t_value * std_error
  
  # Print the results
  cat(mean_value)
  cat("\n95% Confidence Interval: (", ci_lower, ", ", ci_upper, ")\n")
}

ratio <- rep(0,length(beta1))

for (i in 1:length(beta1)){
  mean_value <- mean(smoothed_beta[i,,1,1])
  std_error <- sd(smoothed_beta[i,,1,1])
  t_value <- qnorm(0.975)
  
  # Calculate the 95% confidence interval
  ci_lower <- mean_value - t_value * std_error
  ci_upper <- mean_value + t_value * std_error
  
  if (beta1[i] >= ci_lower & beta1[i] <= ci_upper){
    ratio[i] <- 1
  }
}



# Convert smoothed beta matrix to a data frame for plotting
df_fitted <- data.frame(x1 = x, x2 = z, ind_beta1 = ind_beta1, smo_beta1 = smoothed_beta)

fld <- with(df_fitted, interp(x = x1, y = x2, z = smo_beta1))
fld$x <- c(fld$x, 1.000001)
fld$y <- c(fld$y, 1.000001)
fld$z <- rbind(cbind(fld$z, rep(10, 40)), rep(10,41))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = expression("Gaussian Kernel Smoothed Beta Using " ~ sigma[ker] == 0.045),
               key.title = title(main = "beta1 ", cex.main = 1))

contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(smoothed_beta, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
        main = "Gaussian Kernel Smoothed Beta", xlab = "x1", ylab = "x2")




# III. Gaussian Kernel Smoothing + using original var

# IV. Gaussian Kernel Smoothing + using adjusted var


# Simulation: DD + two betas


# Real Data: HCP Motion






# Example vector of values (replace with your actual data)
data_vector <- rnorm(1000, mean = 5, sd = 2)

# Calculate the sample mean
mean_value <- mean(data_vector)

# Calculate the standard error of the mean
n <- length(data_vector)
std_error <- sd(data_vector) / sqrt(n)

# Get the critical t-value for a 95% confidence level
# Using n-1 degrees of freedom
t_value <- qt(0.975, df = n - 1)

# Calculate the 95% confidence interval
ci_lower <- mean_value - t_value * std_error
ci_upper <- mean_value + t_value * std_error

# Print the results
cat("95% Confidence Interval: [", ci_lower, ", ", ci_upper, "]\n")

