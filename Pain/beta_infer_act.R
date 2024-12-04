# Setting 1:  Activation beta (only 1 and 0) with single parcel



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
library(reshape2)
library(shiny)

# Simulation: DD + one beta


# Model: Canonical DD (two betas)

# I. 2D Square Region

# Define the test function f(x, z)
f_function_1 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (3.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.2)^2 / sigma_x^2) - ((z - 0.3)^2 / sigma_z^2))
  term2 <- (1.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.7)^2 / sigma_x^2) - ((z - 0.8)^2 / sigma_z^2))
  sapply(term1 + term2, function(x) if(x>=3.78){1} else{0})
}

f_function_refer <- function(x, z, sigma_x, sigma_z) {
  term1 <- (3.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.2)^2 / sigma_x^2) - ((z - 0.3)^2 / sigma_z^2))
  term2 <- (1.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.7)^2 / sigma_x^2) - ((z - 0.8)^2 / sigma_z^2))
  term1 + term2
}

# Set parameters
sigma_x <- 0.3
sigma_z <- 0.4
n_points <- 40^2  # number of points to generate

# Generate random (x, z) points within the unit square
set.seed(123)
grid <- expand.grid(x = seq(-0.1, 1.1, length.out = sqrt(n_points)), y = seq(-0.1, 1.1, length.out = sqrt(n_points)))

x <- grid[,1]
z <- grid[,2]

beta1 <- f_function_1(x, z, sigma_x, sigma_z)

# Create a data frame from the vectors
data <- data.frame(x1 = x, x2 = z)

# Create a grid for plotting the true function
grid_size <- 50
x_grid <- seq(-0.1, 1.1, length.out = grid_size)
z_grid <- seq(-0.1, 1.1, length.out = grid_size)
beta1_grid <- outer(x_grid, z_grid, function(x, z) f_function_1(x, z, sigma_x, sigma_z))

# Draw contour plots
# Contour plot of the true function
contour(x_grid, z_grid, beta1_grid, main = "True Beta1", xlab = "x1", ylab = "x2",
        xlim= c(-0.1,1), ylim=c(-0.1,1.1))

beta1_grid_refer <- outer(x_grid, z_grid, function(x, z) f_function_refer(x, z, sigma_x, sigma_z))
contour(x_grid, z_grid, beta1_grid_refer, main = "True Beta1", xlab = "x1", ylab = "x2",
        xlim= c(-0.1,1), ylim=c(-0.1,1), levels = 3.78)

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
       
sim_k <- 10
error_sig <- c(1,2,5,10,20,30)
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
    print(j)
  }
}

save(beta_DD, sig_hat, tc_mat, file = "beta_inf_sim.RData")


# Plot fitted independent beta

for (j in 1:length(error_sig)){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  ind_beta_diff <- ind_beta1 - beta1
  t_value <- (ind_beta1) / apply(sig_hat[,,j,1], 1, mean)
  p_value <- (1 - pt(t_value, sim_k-1))
  t_value_global <- abs((ind_beta1) / mean(apply(sig_hat[,,j,1], 1, mean)))
  df_fitted <- data.frame(x1 = x, x2 = z, ind_beta1 = ind_beta1, t_value = t_value)
  
  #create plot
  
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
  
  contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
          matrix(ind_beta1, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
          main = bquote(hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j])), 
          xlab = "x1", ylab = "x2", levels = 1)
  
  
  contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
          matrix(p_value, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
          main = bquote("p-value of " ~ hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j])),  
          ylab = "x2", levels = c(0.2))
  
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
  
  cat(round(mean_value,2), "  95% Confidence Interval: [", round(ci_lower,2), 
      ", ", round(ci_upper,2), "]\n")
  
}

# SMOOTHING

# II. Gaussian Kernel Smoothing + using smoothing var

# Define bandwidth for the Gaussian kernel
# bandwidth <- diff(x)[1] / (2*sqrt(2*log(2))) #0.02682
bw_vec <- c(0.025, 0.05, 0.075, 0.1)

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
      cat(bw, ", ", j, ", ", k, "\n")
    }
  }
}

save(smoothed_beta, smoothed_Y, file = "beta_inf_sim_smo.RData")


# Inference on sigma tilde hat

mean_value <- mean(sig_hat_tilde)
n <- length(sig_hat_tilde)
std_error <- sd(sig_hat_tilde)

t_value <- qt(0.975, df = n - 1)

ci_lower <- mean_value - t_value * std_error
ci_upper <- mean_value + t_value * std_error

cat(round(mean_value, 2), "95% Confidence Interval: [", round(ci_lower,2), ", ", 
    round(ci_upper,2), "]\n")

# Plot smoothed beta

for (j in 1:length(error_sig)){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  ind_beta_diff <- ind_beta1 - beta1
  t_value <- (ind_beta1) / apply(sig_hat[,,j,1], 1, mean)
  p_value <- (1 - pt(t_value, sim_k-1))

  for (bw in 0:length(bw_vec)){
    png(paste0("/Users/user/Desktop/betafit", j, "-", bw, ".png"), width = 500, height = 400)
    if (bw==0){
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
    png(paste0("/Users/user/Desktop/diff", j, "-", bw, ".png"), width = 500, height = 400)
    if (bw==0){
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
  ind_beta_diff <- ind_beta1 - beta1
  t_value <- (ind_beta1) / apply(sig_hat[,,j,1], 1, mean) * sqrt(sum(X^2))
  p_value <- 2 * (1 - pt(abs(t_value), df = sim_k-1))
  false_pos[j,1,] <- (as.integer(p_value < 0.05) - beta1)
  
  for (bw in 0:length(bw_vec)){
    png(paste0("/Users/user/Desktop/pvalue", j, "-", bw, ".png"), width = 400, height = 400)
    if (bw==0){
      contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
              matrix(p_value, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
              main = bquote("p-value of " ~ hat(beta)[ind] ~ "for" ~ sigma ~ "=" ~ .(error_sig[j])),  
              ylab = "x2", levels = levels, col = colors)
    }else{
      smo_beta1 <- apply(smoothed_beta[,,j,1,bw],1,mean)
      smo_beta1_diff <- smo_beta1 - beta1
      
      for (i in 1:length(beta1)){
        for (k in 1:sim_k){
          sig_hat_tilde[j,i,k,bw] <- sqrt(sum(((smoothed_Y[,i,k,j,bw]-mean(smoothed_Y[,i,k,j,bw])) - 
                                            X * smoothed_beta[i,k,j,1,bw])^2)/(len - 1))
        }
      }
      
      sig_hat_mean <- apply(sig_hat_tilde[j,,,bw],1,mean)
      t_value <- (smo_beta1) / sig_hat_mean * sqrt(sum(X^2))
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
for (j in 1:length(error_sig)){
  for (bw in 0:length(bw_vec)){
    png(paste0("/Users/user/Desktop/FP", j, "-", bw, ".png"), width = 500, height = 400)
    if (bw ==0){
      fp <- false_pos[j,bw+1,]
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
                     main =  bquote("FP and FN for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                      "," ~ h ~ "= 0"),
                     key.title = title(main = "beta1 ", cex.main = 1))
    }else{
      df_fitted <- data.frame(x1 = x, x2 = z, false_pos = false_pos[j,bw+1,])
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
                     main =  bquote("FP and FN for" ~ sigma ~ "=" ~ .(error_sig[j]) ~ 
                                      "," ~ h ~ "=" ~ .(bw_vec[bw])),
                     key.title = title(main = "beta1 ", cex.main = 1))
    }
    dev.off()
  }
}

# FP and FN table
for (i in 1:6){
  for (j in 1:5){
    FP <- sum(false_pos[i,j,]==1) / (length(beta1) - sum(beta1)) * 100
    FN <- sum(false_pos[i,j,]==-1) / sum(beta1) * 100
    cat("FP:", round(FP,2), ", FN:", round(FN,2), "\n")
  }
}

# MSE

for (j in 1:6){
  for (bw in 1:4){
    smo_beta1 <- apply(smoothed_beta[,,j,1,bw],1,mean)
    smo_beta1_diff <- smo_beta1 - beta1
    mse <- mean(smo_beta1_diff^2)
    cat(j, "-", bw, " MSE:", round(mse,4), "\n")
  }
}

for (j in 1:6){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  ind_beta_diff <- ind_beta1 - beta1
  mse <- mean(ind_beta_diff^2)
  cat(j, "-0", " MSE:", round(mse,4), "\n")
}




# TPRS
k_gam <- 50
numbeta <- 1
V <- length(beta1)
X <- Canonical_HRF_fitted$info$X[,-1]
i_perm <- 1:(numbeta * V)
j_perm <- c()
for (i in 1:V){
  j_perm <- c(j_perm, (0:(numbeta-1))*V+i)
}
perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 

tprs_beta <- array(0, dim = c(length(beta1), sim_k, length(error_sig)))
tprs_gamma <- array(0, dim = c(k_gam, sim_k, length(error_sig)))
tprs_sigma <- array(0, dim = c(sim_k, length(error_sig)))

for (j in 1:length(error_sig)){
  ind_beta1 <- apply(beta_DD[,,j,1],1,mean)
  gam_df <- data.frame(beta = ind_beta1, x = x, y = z)
  gamfit <- gam(beta ~ s(x,y, k = k_gam), data = gam_df, method = "REML")
  Blist <- as.matrix(model.matrix(gamfit))
  
  X_ext <- bdiag(replicate(V, X, simplify = FALSE))
  
  B <- bdiag(replicate(numbeta, Blist, simplify = FALSE))
  PB <- perm %*% B
  XtX <- t(X_ext) %*% X_ext
  BPXtX <- t(PB) %*% XtX
  A <- BPXtX %*% PB
  
  for (i in 1:sim_k){
    y_ext <- c(tc_mat[,,i,j])
    Xty_ext <- t(X_ext) %*% y_ext
    b <- t(PB) %*% Xty_ext

    tprs_gamma[,i,j] <- Matrix::solve(A, b)[,1]
    tprs_beta[,i,j] <- (Blist %*% tprs_gamma[1:k_gam,i,j])[,1]
    e <- y_ext - X_ext %*% matrix(tprs_beta[,i,j], ncol=1)
    tprs_sigma[i,j] <- sqrt(sum(e^2) / (length(y_ext) - 1))
  }
}



# Plot tprs beta

for (j in 1:length(error_sig)){
  tprs_beta1 <- apply(tprs_beta[,,j],1,mean)
  tprs_beta_diff <- tprs_beta1 - beta1
  
  png(paste0("/Users/user/Desktop/sim_figure/setting1/tprs/betafit", j, "-50.png"), width = 500, height = 400)
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
                                 "," ~ k[gamma] ~ "= 50"),
                 key.title = title(main = "beta1 ", cex.main = 1))
  dev.off()
  
  png(paste0("/Users/user/Desktop/sim_figure/setting1/tprs/diff", j,"-50.png"), width = 500, height = 400)
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
                                 "," ~ k[gamma] ~ "= 50"),
                 key.title = title(main = "beta1 ", cex.main = 1))
  dev.off()
}

# Plot p-values

X <- Canonical_HRF_fitted$info$X[,2]
false_pos <- array(0,dim = c(length(error_sig), length(beta1)))


for (j in 1:length(error_sig)){
  tprs_beta1 <- apply(tprs_beta[,,j],1,mean)
  t_value <- (tprs_beta1) / mean(tprs_sigma[,j]) * sqrt(sum(X^2))
  p_value <- 2 * (1 - pt(abs(t_value), df = sim_k-1))
  false_pos[j,] <- (as.integer(p_value < 0.05) - beta1)
  
  #png(paste0("/Users/user/Desktop/sim_figure/setting1/tprs/pvalue", j, "-150.png"), width = 400, height = 400)
  #contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
  #        matrix(p_value, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
  #        main = bquote("p-value of " ~ hat(beta)[tprs] ~ "for" ~ sigma ~ 
  #                        "=" ~ .(error_sig[j]) ~ "," ~ k[gamma] ~ "= 150"),  
  #        ylab = "x2", levels = levels, col = colors)
  
  #dev.off()
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

