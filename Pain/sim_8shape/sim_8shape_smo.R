# Access the argument passed from the SLURM script
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

# Setting :  Activation beta (but semi-discrete) with single parcel

# Load Rscripts
# setwd("Pain")
source("../Fit_Canonical_HRF.R")
source("../Fit_Spline.R")
source("../get_parameters2.R")
source("../spm_hrf.R")
source("../Det_Logit.R")
source("../tor_make_deconv_mtx3.R")

# Load necessary libraries

packages <- c( "mgcv", "stats",
               "sp", "tidyverse", "ggplot2",
	       "Matrix", "reshape2")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Simulation: DD + one beta

# I. 2D Square Region

n_points <- 40^2  # number of points to generate

# Generate random (x, z) points within the unit square
set.seed(123)
grid <- expand.grid(x = seq(-0.1, 1.1, length.out = sqrt(n_points)), y = seq(-0.1, 1.1, length.out = sqrt(n_points)))

# Set parameters
sigma_x <- 0.3
sigma_z <- 0.4
n_points <- 40^2  # number of points to generate

x <- grid[,1]
z <- grid[,2]

# Define the test function f(x, z)
f_function_1 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (3.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.25)^2 / sigma_x^2) - ((z - 0.35)^2 / sigma_z^2))
  term2 <- (1.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.75)^2 / sigma_x^2) - ((z - 0.85)^2 / sigma_z^2))
  sapply(term1 + term2, function(x) if(x>=3.78){1} else{0})
}

# Function to compute Gaussian kernel weight
gaussian_weight <- function(dx, dy, bandwidth) {
  exp(- (dx^2 + dy^2) / (2 * bandwidth^2))
}

beta10 <- f_function_1(x, z, sigma_x, sigma_z)

beta1 <- beta10

bandwidth <- 0.08

for (k in 1:length(beta10)){
  gx <- x[k]
  gy <- z[k]
  weights <- mapply(function(x, y) gaussian_weight(gx - x, gy - y, bandwidth), x, z)
  beta1[k] <- sum(weights * beta10) / sum(weights)
}



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

sim_k <- 100
error_sig <- c(1,2,5,8,10,15,20)
#beta_DD <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1))
#sig_hat <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1))
#tc_mat <- array(0, dim = c(len, length(beta1), sim_k, length(error_sig)))

#for (j in 1:length(error_sig)){
#  for (k in 1:length(beta1)){
#    true_sig <-  beta1[k] * conv(Run, CanonicalBasisSet(TR)$h)[1:len]
#    xsecs <- 0:40
   
#    for (i in 1:sim_k){
#      tc_noise <- noise_arp(n = len, phi = c(0, 0), sigma = error_sig[j])
#      tc <- true_sig + tc_noise
#      tc <- true_sig + tc_noise
#      tc_mat[,k,i,j] <- tc
      
      #Canonical HRF
#      Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Runc, T, 1)
#      sig_hat[k,i,j,] <- sqrt(sum(Canonical_HRF_fitted$e^2)/(length(Canonical_HRF_fitted$e)-1))
#      beta_DD[k,i,j,] <- Canonical_HRF_fitted$info$b[1,]
#    }
#  }
#}

#save(beta_DD, sig_hat, tc_mat, file = "beta_inf_sim.RData")
load("beta_inf_sim.RData")

##################################GKS

# Define bandwidth for the Gaussian kernel
# bandwidth <- diff(x)[1] / (2*sqrt(2*log(2))) #0.02682
bw_vec <- seq(0.01, 0.1, length = 10)


# Calculate smoothed beta and Y values

smoothed_beta <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1, length(bw_vec)))
smoothed_Y <- array(0, dim = c(len, length(beta1), sim_k, length(error_sig), length(bw_vec)))
weights_total <- matrix(0, length(beta1), length(beta1))

for (bw in index:index){
  bandwidth <- bw_vec[bw]
  for (j in 1:length(error_sig)){
    for (k in 1:length(beta1)){
      gx <- x[k]
      gy <- z[k]
      # Calculate weights based on the Gaussian kernel
      weights <- mapply(function(x, y) gaussian_weight(gx - x, gy - y, bandwidth), x, z)
      weights_total[k,] <- weights
      
      for (i in 1:sim_k){
        smoothed_beta[k,i,j,1,bw] <- sum(weights * beta_DD[,i,j,1]) / sum(weights)
        for (l in 1:len){
          smoothed_Y[l,k,i,j,bw] <- sum(weights * tc_mat[l,,i,j]) / sum(weights)
        }
      }
      cat(bw, "-", j, "-", k, "\n")
    }
  }
}

save(smoothed_beta, file = paste0("sim_8shape/sim_8shape_smo", index, ".RData"))
