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
library(ggplot2)
library(tidyverse)
library(sp)
library(gstat)
library(raster)
library(stats)
library(mgcv)  # for thin plate regression spline
library(akima)
library(Matrix)
library(reshape2)
library(shiny)

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

beta10 <- f_function_1(x, z, sigma_x, sigma_z)

beta1 <- array(0, dim = c(length(beta10),length(bw_vec_beta)))
bandwidth <- 0.08

for (k in 1:length(beta10)){
  gx <- x[k]
  gy <- z[k]
  weights <- mapply(function(x, y) gaussian_weight(gx - x, gy - y, bandwidth), x, z)
  weights_total[k,] <- weights
  beta1[k] <- sum(weights * beta10) / sum(weights)
}

# Draw contour plots
# Contour plot of the true function
contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(beta1, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T),  
        main = "True beta", ylab = "x2", levels = 0.35)

df_fitted <- data.frame(x1 = x, x2 = z, beta1 = beta1)
fld <- with(df_fitted, interp(x = x1, y = x2, z = beta1))
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2", 
               main =  "True beta",
               key.title = title(main = "beta1 ", cex.main = 1))

# Simulation

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
tc_mg <- c()

sim_k <- 100
error_sig <- c(1,2,5,8,10,15,20)
beta_DD <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1))
sig_hat <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1))
tc_mat <- array(0, dim = c(len, length(beta1), sim_k, length(error_sig)))

for (j in 1:length(error_sig)){
  for (k in 1:length(beta1)){
    true_sig <-  beta1[k] * conv(Run, CanonicalBasisSet(TR)$h)[1:len]
    xsecs <- 0:40
    
    for (i in 1:sim_k){
      tc_noise <- noise_arp(n = len, phi = c(0, 0), sigma = error_sig[j])
      tc <- true_sig + tc_noise
      tc <- true_sig + tc_noise
      tc_mat[,k,i,j] <- tc
      
      #Canonical HRF
      Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Runc, T, 1)
      sig_hat[k,i,j,] <- sqrt(sum(Canonical_HRF_fitted$e^2)/(length(Canonical_HRF_fitted$e)-1))
      beta_DD[k,i,j,] <- Canonical_HRF_fitted$info$b[1,]
    }
  }
}

save(beta_DD, sig_hat, tc_mat, file = "beta_inf_sim.RData")


##################################GKS

# Define bandwidth for the Gaussian kernel
# bandwidth <- diff(x)[1] / (2*sqrt(2*log(2))) #0.02682
bw_vec <- seq(0.01, 0.1, length = 10)

# Function to compute Gaussian kernel weight
gaussian_weight <- function(dx, dy, bandwidth) {
  exp(- (dx^2 + dy^2) / (2 * bandwidth^2))
}

# Calculate smoothed beta and Y values

smoothed_beta <- array(0, dim = c(length(beta1), sim_k, length(error_sig), 1, length(bw_vec)))
smoothed_Y <- array(0, dim = c(len, length(beta1), sim_k, length(error_sig), length(bw_vec)))
weights_total <- matrix(0, length(beta1), length(beta1))

for (bw in 1:length(bw_vec)){
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

save(smoothed_beta, file = "beta_inf_sim_smo.RData")


###################################TPRS

k_gam_vec <- c(20, 40, 60, 80, 100, 120, 140, 160, 180)
numbeta <- 1
V <- length(beta1)
X <- Canonical_HRF_fitted$info$X[,-1]
i_perm <- 1:(numbeta * V)
j_perm <- c()
for (i in 1:V){
  j_perm <- c(j_perm, (0:(numbeta-1))*V+i)
}
X_ext <- bdiag(replicate(V, X, simplify = FALSE))
perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 
tprs_beta <- array(0, dim = c(length(beta1), sim_k, length(error_sig), length(k_gam_vec)))
tprs_gamma <- list() 
array(0, dim = c(k_gam, sim_k, length(error_sig), length(k_gam_vec)))
tprs_sigma <- array(0, dim = c(sim_k, length(error_sig), length(k_gam_vec)))

for (k_ind in 1:length(k_gam_vec)){
  k_gam <- k_gam_vec[k_ind]
  for (j in 1:length(error_sig)){
    ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
    gam_df <- data.frame(beta = ind_beta1, x = x, y = z)
    gamfit <- gam(beta ~ s(x,y, k = k_gam), data = gam_df, method = "REML")
    Blist <- as.matrix(model.matrix(gamfit))
  
    B <- bdiag(replicate(numbeta, Blist, simplify = FALSE))
    PB <- perm %*% B
    XtX <- t(X_ext) %*% X_ext
    BPXtX <- t(PB) %*% XtX
    A <- BPXtX %*% PB
    
    for (i in 1:sim_k){
      y_ext <- c(tc_mat[,,i,j])
      Xty_ext <- t(X_ext) %*% y_ext
      b <- t(PB) %*% Xty_ext
      gam_sol <- Matrix::solve(A, b)
      tprs_gamma[[k_ind]] <- gam_sol[,1]
      tprs_beta[,i,j,k_ind] <- (Blist %*% gam_sol)[,1]
      e <- y_ext - X_ext %*% matrix(tprs_beta[,i,j,k_ind], ncol=1)
      tprs_sigma[i,j,k_ind] <- sqrt(sum(e^2) / (length(y_ext) - 1))
    }
  }
}

save(tprs_beta, file = "beta_inf_sim_tprs.RData")
