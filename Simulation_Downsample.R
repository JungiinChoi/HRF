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

# Read 25 Downsampled HRF
raw_HRF <- read.delim("JunginHRF25_Downsample10.txt", header = FALSE)
HRFdf <- data.frame()
HRFmat <- matrix(nrow = 25, ncol = 640)
for (i in 1:25){
  HRF_tmp <- as.numeric(strsplit(raw_HRF[i,], split = "\\s+")[[1]][-1])
  HRFdf <- rbind(HRFdf, data.frame(index = rep(i,640), x = 0:639, y = HRF_tmp))
  HRFmat[i,] <- HRF_tmp
}

HRFdf$index <- factor(HRFdf$index)

# Plot the hrf
ggplot(HRFdf, aes(x = x/10, y, group = index, colour = index)) +
  xlab("x") +
  geom_line()

# Settings
TR <- 0.5
Down <- 10  # Downsampling factor
T <- 30
b <- 1

# Set onsets at downsampled resolution
set.seed(1001)
R <- sample(1:(640 * Down), 36)
Run <- rep(0, 640 * Down)
Run <- matrix(Run, ncol = 1)
Run[R] <- 1
Runc <- list(Run)
true_sig <- b * conv(Run, hrf)
true_sig <- true_sig[1:(640 * Down)]

onset <- rep(-100,6400)
onset[R] <- R

HRF_df <- data.frame(t = 1:len, tc = tc, onset = onset)

ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc, color = 'Data')) +
  theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 6400)) +
  scale_color_manual(name = "", values=c('blue'))

tc <- tc[seq(1, length(tc), by = Down)] # Put data back into original time resolution

# Settings
numstim <- 1
len <- length(Run)
t <- seq(1, T, by = TR / Down)

K <- 8
norder <- 4

# Downsampling Simulation
sim_list <- list()
for (h in 1:25){
  true_sig <- b * conv(Run, HRFmat[h,])[1:640]
  tc_noise <- noise_arp(n = 640, phi = c(0.3, 0))
  tc <- true_sig + 0.5 * tc_noise
  tc <- (tc - mean(tc)) / sd(tc)
  xsecs <- 0:40
  
  true_sig <- b * conv(Run, HRFmat[h,])
  true_sig <- true_sig[1:(640 * Down)]
  
  tc_mat <- matrix(0, nrow = 640, ncol = 100)

  fitted_hrf_spline <- matrix(0, nrow = 581, ncol = 100)
  params_spline <- matrix(0, nrow = 3, ncol = 100)
  MSE_spline <- matrix(0, nrow = 2, ncol = 100)
  
  for (i in 1:100){
    tc_noise <- noise_arp(n = 640, phi = c(0.3, 0))
    tc <- true_sig + 0.5 * tc_noise
    tc <- (tc - mean(tc)) / sd(tc)
    tc <- tc[seq(1, length(tc), by = Down)]
    xsecs <- 0:580
    tc_mat[,i] <- tc
    
    #B-spline
    Spline_fitted <- Fit_Spline_down(tc, TR, Runc, T, Down)
    fitted_hrf_spline[,i] <- Spline_fitted$hrf[,1]
    params_spline[,i] <- Spline_fitted$param
    e4 <- Spline_fitted$e
    MSE_spline[,i] <- c((1 / (len - 1) * sum(e4^2)), ResidScan(e4, FWHM)$p)
  }
  
  mat_list <- list(tc_mat = tc_mat, fitted_hrf_spline = fitted_hrf_spline, params_spline = params_spline,
                   MSE_spline = MSE_spline)
  sim_list[[h]] <- mat_list
}

# Visualization

## Visualize fitted hrf for each method & each hrf
h = 3
p <- ggplot(data.frame(xsecs = seq(1, 30, by = TR / Down), hrf = HRFmat[h,1:581])) + 
  geom_line(aes(x=xsecs, y = hrf), color = "black", size = 1.5)

mat_list <- sim_list[[h]]

for (i in 1:100){
  hrf_Spline <- mat_list$fitted_hrf_spline[,i]
  p <- p + geom_line(data = data.frame(xsecs = seq(1, 30, by = TR / Down), hrf = hrf_Spline),
                     aes(x = xsecs, y = hrf, color = "Spline"), alpha=0.2)
}
p + scale_color_manual(name = "", values = c("IL" = "red", "sFIR" = "green", "DD" = "magenta", "Spline" = "blue", "NL" = "yellow"))


# Bias

true_params <- t(apply(HRFmat, 1, function(x) get_parameters2(x,1:640))) / 10
true_params[,1] <- rep(1,25)

bias_spline_A <- rep(0,25)

for (h in 1:25){
  mat_list <- sim_list[[h]]
  bias_spline_A[h] <- mean(mat_list$params_spline[3,]/10 - true_params[h,3]) / true_params[h,3]
}
length(bias_spline_A)

IL_A <- data.frame(x = factor(round(true_params[,2])), y = factor(true_params[,3]), 
                   bias = abs(bias_sFIR_A))

ggplot(IL_A, aes(x = x, y = y, fill = bias)) +
  geom_tile(color = "black") +
  coord_fixed() +
  xlab("Time-to-peak") +
  ylab("Width") +
  scale_fill_gradient(low = "white", high = "green", limits = c(0, 2)) 
