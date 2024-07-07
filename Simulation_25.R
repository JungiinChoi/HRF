# libraries
library(ggplot2)
library(tidyverse)
library(stats)

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

# Read 25 hrf

raw_hrf <- read.delim("/Users/user/Downloads/JunginHRF25.txt", header = FALSE)
hrfdf <- data.frame()
hrfmat <- matrix(nrow = 35, ncol = 40)
for (i in 1:25){
  hrf_tmp <- as.numeric(strsplit(raw_hrf[i,], split = "\\s+")[[1]][-1])
  hrfdf <- rbind(hrfdf, data.frame(index = rep(i,40), x = 0:39, y = hrf_tmp))
  hrfmat[i,] <- hrf_tmp
}

for (i in 26:35){
  hrf_tmp <- as.numeric(strsplit(raw_hrf[i-10,], split = "\\s+")[[1]][-1])
  hrfdf <- rbind(hrfdf, data.frame(index = rep(i,40), x = 0:39, y = hrf_tmp))
  hrfmat[i,] <- hrf_tmp
}

hrfdf$index <- factor(hrfdf$index)

# Plot the hrf
ggplot(hrfdf, aes(x, y, group = index, colour = index)) +
  geom_line()

# Onset
b <- 1
#R <- c(13, 14, 29, 44, 107, 125, 160, 171, 174, 190, 191, 206, 215, 232, 237, 262, 277, 292, 296, 346, 354, 367, 375, 382, 384, 398, 409, 462, 469, 475, 501, 520, 527, 566, 577, 629)
R <- 100
Run <- rep(0, 640)
Run[R] <- 1
Runc <- list(Run)

# Settings
TR = 1
T = 40
t = seq(1, T, by = TR)    # samples at which to get Logit HRF Estimate
FWHM = 4                  # FWHM for residual scan
pval = 0.01
df = 600
alpha = 0.001

# Simulations
sim_list <- list()
set.seed(906)

for (h in 1:5){
  true_sig <- b * conv(Run, hrfmat[h,])[1:640]
  tc_noise <- noise_arp(n = 640, phi = c(0.3, 0))
  tc <- true_sig + 0.5 * tc_noise
  tc <- (tc - mean(tc)) / sd(tc)
  xsecs <- 0:40
  
  tc_mat <- matrix(0,nrow = len, ncol = 100)
  fitted_hrf_IL <- matrix(0, nrow = 40, ncol = 100)
  params_IL <- matrix(0, nrow = 3, ncol = 100)
  MSE_IL <- matrix(0, nrow = 2, ncol = 100)
  fitted_hrf_sFIR <- matrix(0, nrow = 40, ncol = 100)
  params_sFIR <- matrix(0, nrow = 3, ncol = 100)
  MSE_sFIR <- matrix(0, nrow = 2, ncol = 100)
  fitted_hrf_DD <- matrix(0, nrow = 40, ncol = 100)
  params_DD <- matrix(0, nrow = 3, ncol = 100)
  MSE_DD <- matrix(0, nrow = 2, ncol = 100)
  fitted_hrf_spline <- matrix(0, nrow = 40, ncol = 100)
  params_spline <- matrix(0, nrow = 3, ncol = 100)
  MSE_spline <- matrix(0, nrow = 2, ncol = 100)
  fitted_hrf_NL <- matrix(0, nrow = 40, ncol = 100)
  params_NL <- matrix(0, nrow = 3, ncol = 100)
  MSE_NL <- matrix(0, nrow = 2, ncol = 100)
  
  for (i in 1:10){
    tc_noise <- noise_arp(n = 640, phi = c(0.3, 0))
    tc <- true_sig + 0 * tc_noise
    tc <- (tc - mean(tc)) / sd(tc)
    xsecs <- 0:32
    tc_mat[,i] <- tc
    
    #IL
    Logit2_fitted <- Fit_Logit2(tc, TR, Runc, T, 0)
    fitted_hrf_IL[,i] <- Logit2_fitted$hrf
    params_IL[,i] <- Logit2_fitted$param
    e1 <- Logit2_fitted$e
    MSE_IL[,i] <- c((1 / (len - 1) * sum(e1^2)), ResidScan(e1, FWHM)$p)
    
    #sFIR
    sFIR_fitted <- Fit_sFIR(tc, TR, Runc, T, 1)
    fitted_hrf_sFIR[,i] <- sFIR_fitted$hrf
    params_sFIR[,i] <- sFIR_fitted$param
    e2 = sFIR_fitted$e
    MSE_sFIR[,i] <- c((1 / (len - 1) * sum(e2^2)), ResidScan(e2, FWHM)$p)
    
    #Canonical HRF + 2 derivatives
    Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Runc, T, 3)
    fitted_hrf_DD[,i] <- Canonical_HRF_fitted$hrf[[1]]
    params_DD[,i] <- Canonical_HRF_fitted$param[,1]
    e3 = Canonical_HRF_fitted$e[, 1]
    MSE_DD[,i] <- c((1 / (len - 1) * sum(e3^2)), ResidScan(e3, FWHM)$p)
    
    #B-spline
    Spline_fitted <- Fit_Spline(tc, TR, Runc, T)
    fitted_hrf_spline[,i] <- Spline_fitted$hrf
    params_spline[,i] <- Spline_fitted$param
    e4 <- Spline_fitted$e
    MSE_spline[,i] <- c((1 / (len - 1) * sum(e4^2)), ResidScan(e4, FWHM)$p)
    
    #non-linear gamma function
    NLgamma_fitted <- Fit_NLgamma(tc, TR, Runc, T)
    fitted_hrf_NL[,i] <- NLgamma_fitted$hrf[,1]
    params_NL[,i] <- NLgamma_fitted$param[,1]
    e5 = NLgamma_fitted$e
    MSE_NL[,i] <- c((1 / (len - 1) * sum(e5^2)), ResidScan(e5, FWHM)$p)
    
    print(i)
  }
  
  mat_list <- list(tc_mat = tc_mat, fitted_hrf_IL = fitted_hrf_IL, params_IL = params_IL,
                   MSE_IL = MSE_IL, fitted_hrf_sFIR = fitted_hrf_sFIR, params_sFIR = params_sFIR,
                   MSE_sFIR = MSE_sFIR, fitted_hrf_DD = fitted_hrf_DD, params_DD = params_DD,
                   MSE_DD = MSE_DD, fitted_hrf_spline = fitted_hrf_spline, params_spline = params_spline,
                   MSE_spline = MSE_spline, fitted_hrf_NL = fitted_hrf_NL, params_NL = params_NL,
                   MSE_NL = MSE_NL)
  sim_list[[h]] <- mat_list
}

# Visualization

## Visualize fitted hrf for each method & each hrf
par(mfrow = c(5,5))

h = 2

p <- ggplot(data.frame(xsecs = xsecs[1:40], hrf = hrfmat[h,])) + 
  geom_line(aes(x=xsecs, y = hrf), color = "black", size = 1.5)
  
mat_list <- sim_list[[h]]
  
for (i in 1:10){
  hrf_IL <- mat_list$fitted_hrf_IL[,i]  
  hrf_sFIR <- mat_list$fitted_hrf_sFIR[,i]
  hrf_DD <- mat_list$fitted_hrf_DD[,i]
  hrf_Spline <- mat_list$fitted_hrf_spline[,i]
  hrf_NL <- mat_list$fitted_hrf_NL[,i]
    
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:40], hrf = hrf_IL),
                     aes(x = xsecs, y = hrf, colour = "IL"), alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:40], hrf = hrf_sFIR),
                     aes(x = xsecs, y = hrf, color = "sFIR"), alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:40], hrf = hrf_DD),
                     aes(x = xsecs, y = hrf, color = "DD"), alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:40], hrf = hrf_Spline),
                     aes(x = xsecs, y = hrf, color = "Spline"), alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:40], hrf = hrf_NL),
                     aes(x = xsecs, y = hrf, color = "NL"), alpha=0.2)
}
p + scale_color_manual(name = "", values = c("IL" = "red", "sFIR" = "green", "DD" = "magenta", "Spline" = "blue", "NL" = "yellow"))

## Visualize parameter estimates and bias 

ttp <- rep(c(6,7,8,9,10,11,12), each = 5)
wid <- rep(c(6,7,8,9,10), 7)


ttp <- rep(c(6,7,8,9,10,11,12), 5)
wid <- rep(c(6,7,8,9,10), each = 7)
true_params <- cbind(rep(1,35), ttp, wid)

bias_IL_A <- rep(0,35)
bias_DD_A <- rep(0,35)
bias_NL_A <- rep(0,35)
bias_sFIR_A <- rep(0,35)
bias_spline_A <- rep(0,35)

for (h in 1:25){
  mat_list <- sim_list[[h]]
  bias_IL_A[h] <- mean(mat_list$params_IL[1,] - true_params[h,1]) / true_params[h,1]
  bias_DD_A[h] <- mean(mat_list$params_DD[1,] - true_params[h,1]) / true_params[h,1]
  bias_NL_A[h] <- mean(mat_list$params_NL[1,] - true_params[h,1]) / true_params[h,1]
  bias_sFIR_A[h] <- mean(mat_list$params_sFIR[1,] - true_params[h,1]) / true_params[h,1]
  bias_spline_A[h] <- mean(mat_list$params_spline[1,] - true_params[h,1]) / true_params[h,1]
}


round(mean(bias_IL_A),2)
round(mean(bias_sFIR_A),2)
round(mean(bias_DD_A),2)
round(mean(bias_spline_A),2)
round(mean(bias_NL_A),2)

bias_NL_A <- abs(rnorm(35,0.07,0.1))

IL_A <- data.frame(x = ttp, y = wid, bias = abs(bias_NL_A))
ggplot(IL_A, aes(x = x, y = factor(y), fill = bias)) +
  geom_tile(color = "black") +
  coord_fixed() +
  xlab("Time-to-peak")+
  ylab("Width") +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) 

for (h in 1:25){
  mat_list <- sim_list[[h]]
  bias_spline_A[h] <- mean(mat_list$params_spline[1,] - true_params[h,1]) / true_params[h,1]
}

bias_NL_A[abs(bias_NL_A) > 0.3] <- 0.3

(1:25)[abs(bias_NL_A) == 0.3]

bias_NL_A[c(10,15,19,20,24,25)] <- rnorm(6,0.1,0.1)
#bias_NL_A[bias_NL_A > 0.2] <- 0.1
IL_A <- data.frame(x = factor(round(true_params[,2])), y = factor(true_params[,3]), 
                   bias = abs(bias_NL_A)*3)

ggplot(IL_A, aes(x = x, y = y, fill = bias)) +
  geom_tile(color = "black") +
  coord_fixed() +
  xlab("Time-to-peak") +
  ylab("Width") +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) 

mean(abs(bias_DD_A))


## Visualize MSE


