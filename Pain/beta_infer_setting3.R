# Setting 3:  Activation beta (but semi-discrete) with single parcel


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

save(beta_DD, sig_hat, tc_mat, file = "beta_inf_sim3.RData")


# Plot fitted independent beta

for (j in 1:length(error_sig)){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  ind_beta_diff <- ind_beta1 - beta1
  df_fitted <- data.frame(x1 = x, x2 = z, ind_beta1 = ind_beta1, ind_beta_diff = ind_beta_diff)
  
  #create plot
  png(paste0("figure/setting3/GKS/betafit", j, "-0.png"), width = 500, height = 400)
  fld <- with(df_fitted, interp(x = x1, y = x2, z = ind_beta1))
  filled.contour(x = fld$x,
                 y = fld$y,
                 z = fld$z,
                 color.palette =
                   colorRampPalette(c("white", "blue")),
                 xlab = "X1",
                 ylab = "X2",
                 main = bquote(hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j])),
                 key.title = title(main = "beta1 ", cex.main = 1))
  dev.off()
  png(paste0("figure/setting3/GKS/diff", j, "-0.png"), width = 500, height = 400)
  fld <- with(df_fitted, interp(x = x1, y = x2, z = ind_beta_diff))
  filled.contour(x = fld$x,
                 y = fld$y,
                 z = fld$z,
                 color.palette =
                   colorRampPalette(c("red", "white", "blue")),
                 xlab = "X1",
                 ylab = "X2",
                 main = bquote(hat(beta)[ind] ~ "-" ~ beta ~ "for" ~ sigma ~ "=" ~ .(error_sig[j])),
                 key.title = title(main = "beta1 ", cex.main = 1))
  dev.off()
}

# Inference

# Inference on sigma_k

for (j in 1:length(error_sig)){
  mean_value <- mean(apply(sig_hat[,,j,1], 1, mean))
  
  n <- length(apply(sig_hat[,,j,1], 1, mean))
  std_error <- sd(apply(sig_hat[,,j,1], 1, mean))
  
  t_value <- qt(0.975, df = n - 1)
  
  ci_lower <- mean_value - t_value * std_error
  ci_upper <- mean_value + t_value * std_error
  
  #cat(round(mean_value,2), "  95% Confidence Interval: [", round(ci_lower,2), 
  #    ", ", round(ci_upper,2), "]\n")
  
}


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

for (bw in 9:10){
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

save(smoothed_beta, file = "beta_inf_sim3_smo.RData")


# Inference on sigma tilde hat

mean_value <- mean(sig_hat_tilde)
n <- length(sig_hat_tilde)
std_error <- sd(sig_hat_tilde)

t_value <- qt(0.975, df = n - 1)

ci_lower <- mean_value - t_value * std_error
ci_upper <- mean_value + t_value * std_error

#cat(round(mean_value, 2), "95% Confidence Interval: [", round(ci_lower,2), ", ", 
#    round(ci_upper,2), "]\n")

# Plot smoothed beta

for (j in 1:length(error_sig)){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  ind_beta_diff <- ind_beta1 - beta1
  t_value <- (ind_beta1) / apply(sig_hat[,,j,1], 1, mean)
  p_value <- (1 - pt(t_value, sim_k-1))
  
  for (bw in 0:length(bw_vec)){
    if (bw==0){
      png(paste0("figure/setting3/GKS/betafit", j, "-0.00.png"), width = 500, height = 400)
      df_fitted <- data.frame(x1 = x, x2 = z, ind_beta1 = ind_beta1)
      fld <- with(df_fitted, interp(x = x1, y = x2, z = ind_beta1))
      filled.contour(x = fld$x,
                     y = fld$y,
                     z = fld$z,
                     color.palette =
                       colorRampPalette(c("white", "blue")),
                     xlab = "X1",
                     ylab = "X2",
                     main = bquote(hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                     "," ~ h ~ "= 0"),
                     key.title = title(main = "beta1 ", cex.main = 1))
    }else{
      if (bw == 10){
        png(paste0("figure/setting3/GKS/betafit", j, "-0.10.png"), width = 500, height = 400)
      }else{
        png(paste0("figure/setting3/GKS/betafit", j, "-", bw_vec[bw], ".png"), width = 500, height = 400)
      }
      smo_beta1 <- apply(smoothed_beta[,,j,1,bw],1,mean)
      smo_beta1_diff <- smo_beta1 - beta1
      
      df_fitted <- data.frame(x1 = x, x2 = z, smo_beta1 = smo_beta1)
      fld <- with(df_fitted, interp(x = x1, y = x2, z = smo_beta1))
      filled.contour(x = fld$x,
                     y = fld$y,
                     z = fld$z,
                     color.palette =
                       colorRampPalette(c("white", "blue")),
                     xlab = "X1",
                     ylab = "X2",
                     main =  bquote(hat(tilde(beta)) ~ "for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                      "," ~ h ~ "=" ~ .(bw_vec[bw])),
                     key.title = title(main = "beta1 ", cex.main = 1))
    }
    dev.off()
  }
  
  for (bw in 0:length(bw_vec)){
    if (bw==0){
      png(paste0("figure/setting3/GKS/diff", j, "-0.00.png"), width = 500, height = 400)
      df_fitted <- data.frame(x1 = x, x2 = z, ind_beta_diff = ind_beta_diff)
      fld <- with(df_fitted, interp(x = x1, y = x2, z = ind_beta_diff))
      filled.contour(x = fld$x,
                     y = fld$y,
                     z = fld$z,
                     color.palette =
                       colorRampPalette(c("red", "white", "blue")),
                     xlab = "X1",
                     ylab = "X2",
                     main = bquote(hat(beta)[ind] - beta ~ "for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                     "," ~ h ~ "= 0"),
                     key.title = title(main = "beta1 ", cex.main = 1))
    }else{
      if (bw == 10){
        png(paste0("figure/setting3/GKS/diff", j, "-0.10.png"), width = 500, height = 400)
      }else{
        png(paste0("figure/setting3/GKS/diff", j, "-", bw_vec[bw], ".png"), width = 500, height = 400)
      }
      smo_beta1 <- apply(smoothed_beta[,,j,1,bw],1,mean)
      smo_beta1_diff <- smo_beta1 - beta1
      
      df_fitted <- data.frame(x1 = x, x2 = z, smo_beta1_diff = smo_beta1_diff)
      fld <- with(df_fitted, interp(x = x1, y = x2, z = smo_beta1_diff))
      filled.contour(x = fld$x,
                     y = fld$y,
                     z = fld$z,
                     color.palette =
                       colorRampPalette(c("red", "white", "blue")),
                     xlab = "X1",
                     ylab = "X2",
                     main =  bquote(hat(tilde(beta)) - beta ~ "for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                      "," ~ h ~ "=" ~ .(bw_vec[bw])),
                     key.title = title(main = "beta1 ", cex.main = 1))
    }
    dev.off()
  }
}

# Plot p-values

sig_hat_tilde <- array(0, dim = c(length(error_sig), length(beta1), sim_k, length(bw_vec)))
X <- Canonical_HRF_fitted$info$X[,2]
levels <- c(0.01,0.03,0.05,0.08,0.1,0.15,0.2,0.3)
colors <- c(rep("black",2), "red", rep("black",5))
false_pos <- array(0,dim = c(length(error_sig), length(bw_vec)+1, length(beta1)))

for (j in 1:length(error_sig)){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  ind_beta_sig <- apply(beta_DD[,,j,1],1,sd)
  ind_beta_diff <- ind_beta1 - beta1
  t_value <- (ind_beta1) / ind_beta_sig
  p_value <- 2 * (1 - pt(abs(t_value), df = sim_k-1))
  false_pos[j,1,] <- (as.integer(p_value < 0.05) - beta1)
  
  for (bw in 0:length(bw_vec)){
    if (bw==0){
      png(paste0("figure/setting3/GKS/pvalue", j, "-0.00.png"), width = 400, height = 400)
      contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
              matrix(p_value, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
              main = bquote("p-value of " ~ hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j])),  
              ylab = "x2", levels = levels, col = colors)
    }else{
      if (bw == 10){
        png(paste0("figure/setting3/GKS/pvalue", j, "-0.10.png"), width = 400, height = 400)
      }else{
        png(paste0("figure/setting3/GKS/pvalue", j, "-", bw_vec[bw], ".png"), width = 400, height = 400)
      }
      smo_beta1 <- apply(smoothed_beta[,,j,1,bw],1,mean)
      smo_beta1_sig <- apply(smoothed_beta[,,j,1,bw],1,sd)
      smo_beta1_diff <- smo_beta1 - beta1
      
      for (i in 1:length(beta1)){
        for (k in 1:sim_k){
          sig_hat_tilde[j,i,k,bw] <- sqrt(sum(((smoothed_Y[,i,k,j,bw]-mean(smoothed_Y[,i,k,j,bw])) - 
                                                 X * smoothed_beta[i,k,j,1,bw])^2)/(len - 1))
        }
      }

      t_value <- (smo_beta1) / smo_beta1_sig
      p_value <- 2 * (1 - pt(abs(t_value), df = sim_k-1))
      false_pos[j,bw+1,] <- (as.integer(p_value < 0.05) - beta1)
      contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
              matrix(p_value, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
              main = bquote("p-value of " ~ hat(tilde(beta)) ~ "for" ~ sigma ~ 
                              "=" ~ .(error_sig[j]) ~ "," ~ h ~ "=" ~ .(bw_vec[bw])),  
              ylab = "x2", levels = levels, col = colors)
    }
    dev.off()
  }
}


# Plot FP(False positive) and FN(False negative)

ind_gks_beta <- array(beta_DD[,,,1], dim = c(1600,10,7,1))
gks_beta <- abind::abind(ind_gks_beta, smoothed_beta[,,,1,], along = 4) # 1600 10 7 11
beta_10 <- sapply(beta1, function(x)if (x>=0.35){1}else{0})
bw_vec_tmp <- c(0,bw_vec)
fpfn_array <- array(0, dim = c(length(beta1), length(error_sig), length(bw_vec_tmp)))
for (j in 1:length(error_sig)){
  for (bw in 0:length(bw_vec)){
    if (bw == 0){
      png(paste0("figure/setting3/GKS/FPFN", j, "-0.00.png"), width = 500, height = 400)
    }else if (bw == 10){
      png(paste0("figure/setting3/GKS/FPFN", j, "-0.10.png"), width = 500, height = 400)
    }else{
      png(paste0("figure/setting3/GKS/FPFN", j, "-", bw_vec[bw], ".png"), width = 500, height = 400)
    }
    beta_mean_tmp <- apply(gks_beta[,,j,bw+1],1,mean)
    beta_sd_tmp <- apply(gks_beta[,,j,bw+1],1,sd)
    t_value <- (beta_mean_tmp) / beta_sd_tmp
    p_value <- 2 * (1 - pt(abs(t_value), df = sim_k-1))
    fpfn <- sapply(p_value, function(x) if (x<=0.05){1}else{0}) - beta_10
    fpfn_array[,j,bw+1] <- fpfn 
    df_fitted <- data.frame(x1 = x, x2 = z, fpfn = fpfn)
    fld <- with(df_fitted, interp(x = x1, y = x2, z = fpfn))
    filled.contour(x = fld$x,
                   y = fld$y,
                   z = fld$z,
                   color.palette =
                     colorRampPalette(c("yellow", "white", "red")),
                   levels = c(-1.1, -0.9, -0.1, 0.1, 0.9, 1.1),
                   col = c("yellow","yellow","white", "red", "red", "red"),
                   xlab = "X1",
                   ylab = "X2",
                   main =  bquote("FP and FN for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                    "," ~ h ~ "=" ~ .(bw_vec_tmp[bw+1])),
                   key.title = title(main = "beta1 ", cex.main = 1))
    dev.off()
  }
}


fpfn_table <- data.frame(sigma = c(), bw = c(), fp = c(), fn = c())
for (i in 1:length(error_sig)){
  for (j in 1:length(bw_vec_tmp)){
    fp_tmp <- sum(fpfn_array[,i,j] == 1) / sum(beta_10 == 0) * 100
    fn_tmp <- sum(fpfn_array[,i,j] == -1) / sum(beta_10 == 1) * 100
    fpfn_table <- rbind(fpfn_table, data.frame(sigma = error_sig[i], bw = bw_vec_tmp[j],
                                               fp = fp_tmp, fn = fn_tmp))
  }
}



# MSE

for (j in 1:7){
  for (bw in 1:10){
    smo_beta1 <- apply(smoothed_beta[,,j,1,bw],1,mean)
    smo_beta1_diff <- smo_beta1 - beta1
    mse <- mean(smo_beta1_diff^2)
    cat(j, "-", bw, " MSE:", round(mse*100,2), "\n")
  }
}

for (j in 1:7){
  for (k_gam in 1:length(k_gam_vec)){
    tprs_beta1_tmp <- apply(tprs_beta[,,j,k_gam],1,mean)
    tprs_beta1_diff_tmp <- tprs_beta1_tmp - beta1
    mse <- mean(tprs_beta1_diff_tmp^2)
    cat(j, "-", k_gam, " MSE:", round(mse*100,2), "\n")
  }
}

for (j in 1:7){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  ind_beta_diff <- ind_beta1 - beta1
  mse <- mean(ind_beta_diff^2)
  cat(j, "-0", " MSE:", round(mse*100,2), "\n")
}

# Optimal Bandwidth Visualization

## Optimal bandwidth graph for specific voxel

data_matrix <- diff_gks_beta[530,1,] # 50 82 88 96 101 1107 1113 1114
bias <- (mean_gks_beta[530,1,] - exp_smoothed_beta[530,])^2
variance <- data_matrix - bias
plot(x_values, data_matrix, type = "b", pch = 19, col = "blue",
     xlab = "Gaussian Kernel Bandwidth", ylab = "", 
     main = bquote("||"~hat(tilde(beta)) - beta~"||"^2~"for (0.18, 0.30), " ~ sigma == .(error_sig[1])),
     ylim = c(0,0.006))

lines(x_values, bias, type = "b", col = "red", pch = 19) # Add second vector
lines(x_values, variance, type = "b", col = "green", pch = 19) # Add third vector

# Add legend
legend("topleft", legend = c("MSE", "Bias^2", "Variance"),
       col = c("blue", "red", "green"), lty = 1)


## Optimal bandwidth for each voxel
mean_smo_beta <- apply(smoothed_beta[,,,1,], c(1,3,4), mean) # 1600 7 10
mean_ind_beta <- array(apply(beta_DD[,,,1], c(1,3), mean), dim = c(1600, 7, 1)) # 1600 7
mean_gks_beta <- abind::abind(mean_ind_beta, mean_smo_beta, along = 3) # 1600 7 11

diff_gks_beta <- (mean_gks_beta - array(beta1, dim = c(1600, 7, 11)))^2
optimal_bw <- (apply(diff_gks_beta, 1:2, which.min) - 1) / 100

for (j in 1:length(error_sig)){
  png(paste0("figure/setting3/GKS/opt", j, "-.png"), width = 500, height = 400)
  df_fitted <- data.frame(x1 = x, x2 = z, optimal_bw = optimal_bw[,j])
  fld <- with(df_fitted, interp(x = x1, y = x2, z = optimal_bw))
  filled.contour(x = fld$x,
                 y = fld$y,
                 z = fld$z,
                 color.palette =
                   colorRampPalette(c("white", "blue")),
                 xlab = "X1",
                 ylab = "X2",
                 main = bquote("Optimal kernel bandwidth for" ~ sigma ~ "=" ~ .(error_sig[j])),
                 key.title = title(main = "beta1 ", cex.main = 1))
  dev.off()
}

## Boxplot

colnames(optimal_bw) <- c(1, 2, 5, 8, 10, 15, 20)

optimal_bandwidth_df <- as.data.frame(optimal_bw)

optimal_bandwidth_long <- reshape2::melt(optimal_bandwidth_df, 
                                         variable.name = "sigma", 
                                         value.name = "OptimalBandwidth")

optimal_bandwidth_long$sigma <- as.numeric(as.character(optimal_bandwidth_long$sigma))

ggplot(optimal_bandwidth_long, aes(x = factor(sigma), y = OptimalBandwidth)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.color = "red") +
  labs(
    title = "Boxplot of Optimal Bandwidths for Different Sigma Values",
    x = expression(sigma), 
    y = "Optimal Bandwidth"
  ) +
  theme_minimal()


## Overall optimal bandwidth among 1600 voxels
x_values <- seq(0, 0.1, length = 11)
for (j in 1:length(error_sig)){
  png(paste0("figure/setting3/GKS/overall_opt", j, "-.png"), width = 400, height = 400)
  data_matrix <- diff_gks_beta[,j,]
  mean_values <- colMeans(data_matrix)
  std_error <- apply(data_matrix, 2, sd) / sqrt(nrow(data_matrix))
  ci_upper <- mean_values + 1.96 * std_error  
  ci_lower <- mean_values - 1.96 * std_error  
  
  plot(x_values, mean_values, type = "b", pch = 19, col = "blue", ylim = range(c(ci_lower, ci_upper)),
       xlab = "Gaussian Kernel Bandwidth", ylab = "", 
       main = bquote("||"~hat(tilde(beta)) - beta~"||"^2~"for"~sigma == .(error_sig[j])))
  polygon(c(x_values, rev(x_values)), c(ci_upper, rev(ci_lower)), col = rgb(0, 0, 1, 0.2), border = NA)
  arrows(x_values, ci_lower, x_values, ci_upper, length = 0.05, angle = 90, code = 3, col = "blue")
  dev.off()
  
}

## Cohen's d

compute_cohen_d <- function(beta_estimates) {
  mean_beta <- mean(beta_estimates)
  sd_beta <- sd(beta_estimates)
  
  if (sd_beta > 0) {
    return(mean_beta / sd_beta)
  } else {
    return(NA)
  }
}

cohen_d <- apply(diff_gks_beta_all, c(1, 3, 4), function(slice) compute_cohen_d(slice))

classify_values <- function(x) {
  cut(
    x,
    breaks = c(-Inf, -0.8, -0.5, 0.5, 0.8, Inf),
    labels = c(-2, -1, 0, 1, 2),
    right = FALSE
  )
}

class <- array(as.numeric(classify_values(as.vector(cohen_d))), dim = dim(cohen_d)) - 3

bw_vec_tmp <- c(0,bw_vec)
for (j in 1:length(error_sig)){
  for (bw in 0:length(bw_vec)){
    if (bw==0){
      png(paste0("figure/setting3/GKS/cohen", j, "-0.00.png"), width = 500, height = 400)
    }else{
      if (bw == 10){
        png(paste0("figure/setting3/GKS/cohen", j, "-0.10.png"), width = 500, height = 400)
      }else{
        png(paste0("figure/setting3/GKS/cohen", j, "-", bw_vec[bw], ".png"), width = 500, height = 400)
      }
    }
    df_fitted <- data.frame(x1 = x, x2 = z, class = class[,j,bw+1])
    fld <- with(df_fitted, interp(x = x1, y = x2, z = class))
    filled.contour(x = fld$x,
                   y = fld$y,
                   z = fld$z,
                   color.palette =
                     colorRampPalette(c("yellow","white", "red")),
                   xlab = "X1",
                   ylab = "X2",
                   main =  bquote("Classification by Cohen's d for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                    "," ~ h ~ "=" ~ .(bw_vec_tmp[bw+1])),
                   key.title = title(main = "class", cex.main = 1))
    dev.off()
  }
}

significant_class_alt <- ifelse(abs(beta1) > 0.2, 1, 0)


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

save(tprs_beta, file = "beta_inf_sim3_tprs.RData")

tprs_bic <- log(tprs_sigma^2 * (length(y_ext) - 1) / length(y_ext)) * length(y_ext) 
tprs_bic_mean <- apply(tprs_bic, 2:3, mean)

for (i in 1:7){
  tprs_bic_mean_tmp <- tprs_bic_mean 
  cat(which.min(tprs_bic_mean_tmp[i,]), "\n")
}


# Plot tprs beta
for (k_ind in 1:length(k_gam_vec)){
  for (j in 1:length(error_sig)){
    tprs_beta1 <- apply(tprs_beta[,,j,k_ind],1,mean)
    tprs_beta_diff <- tprs_beta1 - beta1
    png(paste0("figure/setting3/TPRS/betafit", j, "-", k_ind, ".png"), width = 500, height = 400)
    df_fitted <- data.frame(x1 = x, x2 = z, tprs_beta1 = tprs_beta1)
    fld <- with(df_fitted, interp(x = x1, y = x2, z = tprs_beta1))
    filled.contour(x = fld$x,
                   y = fld$y,
                   z = fld$z,
                   color.palette =
                     colorRampPalette(c("white", "blue")),
                   xlab = "X1",
                   ylab = "X2",
                   main = bquote(hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                   "," ~ k[gamma] ~ "=" ~  .(k_gam_vec[k_ind])),
                   key.title = title(main = "beta1 ", cex.main = 1))
    dev.off()
    
    png(paste0("figure/setting3/TPRS/diff", j, "-", k_ind, ".png"), width = 500, height = 400)
    df_fitted <- data.frame(x1 = x, x2 = z, tprs_beta_diff = tprs_beta_diff)
    fld <- with(df_fitted, interp(x = x1, y = x2, z = tprs_beta_diff))
    filled.contour(x = fld$x,
                   y = fld$y,
                   z = fld$z,
                   color.palette =
                     colorRampPalette(c("red","white", "blue")),
                   xlab = "X1",
                   ylab = "X2",
                   main = bquote(hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                   "," ~ k[gamma] ~ "=" ~  .(k_gam_vec[k_ind])),
                   key.title = title(main = "beta1 ", cex.main = 1))
    dev.off()
  }
}


# Plot p-values

X <- Canonical_HRF_fitted$info$X[,2]
false_pos <- array(0,dim = c(length(error_sig), length(beta1)))
levels <- c(0.01,0.03,0.05,0.08,0.1,0.15,0.2,0.3)
colors <- c(rep("black",2), "red", rep("black",5))

for (k_ind in 1:length(k_gam_vec)){
  for (j in 1:length(error_sig)){
    tprs_beta1 <- apply(tprs_beta[,,j,k_ind],1,mean)
    tprs_beta1_sig <- apply(tprs_beta[,,j,k_ind],1,sd)
    t_value <- (tprs_beta1) / tprs_beta1_sig
    p_value <- 2 * (1 - pt(abs(t_value), df = sim_k-1))
    false_pos[j,] <- (as.integer(p_value < 0.05) - beta1)
    
    png(paste0("figure/setting3/TPRS/pvalue", j, "-", k_ind, ".png"), width = 400, height = 400)
    contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
            matrix(p_value, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
            main = bquote("p-value of " ~ hat(beta)[tprs] ~ "for" ~ sigma ~ 
                            "=" ~ .(error_sig[j]) ~ "," ~ k[gamma] ~ "=" ~  .(k_gam_vec[k_ind])),  
            ylab = "x2", levels = levels, col = colors)
    dev.off()
  }
}



# Plot FP(False positive) and FN(False negative)
for (j in 1:length(error_sig)){
  png(paste0("/Users/user/Desktop/sim_figure/setting1/tprs/FP", j, "-50.png"), width = 500, height = 400)
  fp <- false_pos[j,]
  df_fitted <- data.frame(x1 = x, x2 = z, false_pos = fp)
  fld <- with(df_fitted, interp(x = x1, y = x2, z = false_pos))
  filled.contour(x = fld$x,
                 y = fld$y,
                 z = fld$z,
                 color.palette =
                   colorRampPalette(c("yellow", "white", "red")),
                 levels = c(-1.1, -0.9, -0.1, 0.1, 0.9, 1.1),
                 col = c("yellow","white","white", "white", "red", "red"),
                 xlab = "X1",
                 ylab = "X2",
                 main =  bquote("FP and FN for" ~ sigma ~ "=" ~ .(error_sig[j])
                                ~ "," ~ k[gamma] ~ "=50"),
                 key.title = title(main = "beta1 ", cex.main = 1))
  
  dev.off()
}

# FP and FN table
for (i in 1:6){
  FP <- sum(false_pos[i,]==1) / (length(beta1) - sum(beta1)) * 100
  FN <- sum(false_pos[i,]==-1) / sum(beta1) * 100
  cat(k_gam, " FP:", round(FP,2), ", FN:", round(FN,2), "\n")
}

# MSE
for (j in 1:6){
  tprs_beta1 <- apply(tprs_beta[,,j],1,mean)
  tprs_beta_diff <- tprs_beta1 - beta1
  mse <- mean(tprs_beta_diff^2)
  cat(k_gam, " MSE:", round(mse,4), "\n")
}


# Optimal number of basis Visualization

mean_tprs_beta <- apply(tprs_beta[,,,], c(1,3,4), mean) # 1600 7 9
diff_tprs_beta <- (mean_tprs_beta - array(beta1, dim = c(1600, 7, 9)))^2
optimal_kgam <- apply(diff_tprs_beta, 1:2, function(x) k_gam_vec[which.min(x)])

## Optimal number of basis for each voxel

for (j in 1:length(error_sig)){
  png(paste0("figure/setting3/TPRS/opt", j, "-.png"), width = 500, height = 400)
  df_fitted <- data.frame(x1 = x, x2 = z, optimal_kgam = optimal_kgam[,j])
  fld <- with(df_fitted, interp(x = x1, y = x2, z = optimal_kgam))
  filled.contour(x = fld$x,
                 y = fld$y,
                 z = fld$z,
                 color.palette =
                   colorRampPalette(c("white", "blue")),
                 xlab = "X1",
                 ylab = "X2",
                 main = bquote("Optimal number of basis for" ~ sigma ~ "=" ~ .(error_sig[j])),
                 key.title = title(main = "beta1 ", cex.main = 1))
  dev.off()
}

## Boxplot

# Assign sigma values as column names
colnames(optimal_kgam) <- c(1, 2, 5, 8, 10, 15, 20)

# Convert the matrix into a data frame suitable for ggplot2
optimal_kgam_df <- as.data.frame(optimal_kgam)

# Reshape the data into long format
optimal_bandwidth_long <- reshape2::melt(optimal_kgam_df, 
                                         variable.name = "sigma", 
                                         value.name = "OptimalKgam")

# Convert sigma column to numeric
optimal_bandwidth_long$sigma <- as.numeric(as.character(optimal_bandwidth_long$sigma))

# Create the boxplot
ggplot(optimal_bandwidth_long, aes(x = factor(sigma), y = OptimalKgam)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.color = "red") +
  labs(
    title = "Boxplot of Optimal Number of Basis for Different Sigma Values",
    x = expression(sigma), 
    y = "Optimal Number of Basis"
  ) +
  theme_minimal()


## Overall optimal bandwidth among 1600 voxels
x_values <- k_gam_vec
for (j in 1:length(error_sig)){
  png(paste0("figure/setting3/TPRS/overall_opt", j, "-.png"), width = 400, height = 400)
  data_matrix <- diff_tprs_beta[,j,]
  mean_values <- colMeans(data_matrix)
  std_error <- apply(data_matrix, 2, sd) / sqrt(nrow(data_matrix))
  ci_upper <- mean_values + 1.96 * std_error  
  ci_lower <- mean_values - 1.96 * std_error  
  
  plot(x_values, mean_values, type = "b", pch = 19, col = "blue", ylim = range(c(ci_lower, ci_upper)),
       xlab = "Number of Thin Plate Spline Basis", ylab = "", 
       main = bquote("||"~hat(tilde(beta)) - beta~"||"^2~"for"~sigma == .(error_sig[j])))
  polygon(c(x_values, rev(x_values)), c(ci_upper, rev(ci_lower)), col = rgb(0, 0, 1, 0.2), border = NA)
  arrows(x_values, ci_lower, x_values, ci_upper, length = 0.05, angle = 90, code = 3, col = "blue")
  dev.off()
  
}

# Plot FP(False positive) and FN(False negative)

beta_10 <- sapply(beta1, function(x)if (x>=0.35){1}else{0})
fpfn_array <- array(0, dim = c(length(beta1), length(error_sig), length(k_gam_vec)))

for (j in 1:length(error_sig)){
  for (k_gam in 1:length(k_gam_vec)){
    png(paste0("figure/setting3/TPRS/FPFN", j, "-", k_gam, ".png"), width = 500, height = 400)
    beta_mean_tmp <- apply(tprs_beta[,,j,k_gam],1,mean)
    beta_sd_tmp <- apply(tprs_beta[,,j,k_gam],1,sd)
    t_value <- (beta_mean_tmp) / beta_sd_tmp
    p_value <- 2 * (1 - pt(abs(t_value), df = sim_k-1))
    fpfn <- sapply(p_value, function(x) if (x<=0.05){1}else{0}) - beta_10
    fpfn_array[,j,k_gam] <- fpfn 
    df_fitted <- data.frame(x1 = x, x2 = z, fpfn = fpfn)
    fld <- with(df_fitted, interp(x = x1, y = x2, z = fpfn))
    filled.contour(x = fld$x,
                   y = fld$y,
                   z = fld$z,
                   color.palette =
                     colorRampPalette(c("yellow", "white", "red")),
                   levels = c(-1.1, -0.9, -0.1, 0.1, 0.9, 1.1),
                   col = c("yellow","yellow","white", "red", "red", "red"),
                   xlab = "X1",
                   ylab = "X2",
                   main =  bquote("FP and FN for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                    "," ~ k[gamma] ~ "=" ~ .(k_gam_vec[k_gam])),
                   key.title = title(main = "beta1 ", cex.main = 1))
    dev.off()
  }
}

fpfn_table <- data.frame(sigma = c(), k_gam = c(), fp = c(), fn = c())
for (i in 1:length(error_sig)){
  for (j in 1:length(k_gam_vec)){
    fp_tmp <- sum(fpfn_array[,i,j] == 1) / sum(beta_10 == 0) * 100
    fn_tmp <- sum(fpfn_array[,i,j] == -1) / sum(beta_10 == 1) * 100
    fpfn_table <- rbind(fpfn_table, data.frame(sigma = error_sig[i], k_gam = k_gam_vec[j],
                                               fp = fp_tmp, fn = fn_tmp))
  }
}

## Cohen's d

compute_cohen_d <- function(beta_estimates) {
  mean_beta <- mean(beta_estimates)
  sd_beta <- sd(beta_estimates)
  
  if (sd_beta > 0) {
    return(mean_beta / sd_beta)
  } else {
    return(NA)
  }
}

tprs_beta_diff <- tprs_beta - array(beta1, dim = dim(tprs_beta))
cohen_d_tprs <- apply(tprs_beta_diff, c(1, 3, 4), function(slice) compute_cohen_d(slice))

classify_values <- function(x) {
  cut(
    x,
    breaks = c(-Inf, -0.8, -0.5, 0.5, 0.8, Inf),
    labels = c(-2, -1, 0, 1, 2),
    right = FALSE
  )
}

class_tprs <- array(as.numeric(classify_values(as.vector(cohen_d_tprs))), dim = dim(cohen_d_tprs)) - 3

for (j in 1:length(error_sig)){
  for (k_gam in 1:length(k_gam_vec)){
    png(paste0("figure/setting3/TPRS/cohen", j, "-", k_gam, ".png"), width = 500, height = 400)
    
    df_fitted <- data.frame(x1 = x, x2 = z, class = class_tprs[,j,k_gam])
    fld <- with(df_fitted, interp(x = x1, y = x2, z = class))
    filled.contour(x = fld$x,
                   y = fld$y,
                   z = fld$z,
                   color.palette =
                     colorRampPalette(c("yellow","white", "red")),
                   xlab = "X1",
                   ylab = "X2",
                   main =  bquote("Classification by Cohen's d for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                    "," ~ k[gamma] ~ "=" ~ .(k_gam_vec[k_gam])),
                   key.title = title(main = "class", cex.main = 1))
    dev.off()
  }
}

#################### Overall Comparison

exp_smoothed_beta <- array(0, dim = c(1600, 11))
exp_smoothed_beta[,1] <- beta1

for (bw in 1:length(bw_vec)){
  bandwidth <- bw_vec[bw]
  for (k in 1:length(beta1)){
    gx <- x[k]
    gy <- z[k]
    # Calculate weights based on the Gaussian kernel
    weights <- mapply(function(x, y) gaussian_weight(gx - x, gy - y, bandwidth), x, z)
    weights_total[k,] <- weights
    exp_smoothed_beta[k,bw+1] <- sum(weights * beta1) / sum(weights)
    cat(bw, "-", k, "\n")
  }
}

ind_tmp <- which(optimal_bw[,1]==0)
plot(x[ind_tmp], z[ind_tmp], xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "X", ylab = "Y", 
     main = bquote("Optimal Bandwidth = 0 for " ~ sigma ~ "=" ~ .(error_sig[1])),
     pch = 19, col = "blue")

# Optionally add a grid for better visualization
grid(col = "gray", lty = "dotted")


