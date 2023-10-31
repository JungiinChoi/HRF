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
ggplot(HRFdf, aes(x, y, group = index, colour = index)) +
  geom_line()

# Settings
TR <- 0.5
Down <- 10  # Downsampling factor
T <- 30
b <- 1

# Load downsampled HRF
hrf <- HRFmat[2,]

# Set onsets at downsampled resolution
set.seed(1001)
R <- sample(1:(640 * Down), 36)
Run <- rep(0, 640 * Down)
Run <- matrix(Run, ncol = 1)
Run[R] <- 1
Runc <- list(Run)
true_sig <- b * conv(Run, hrf)
true_sig <- true_sig[1:(640 * Down)]

tc_noise <- noise_arp(640 * Down, c(0.3, 0))
tc <- true_sig + 0.5 * tc_noise

tc <- tc[seq(1, length(tc), by = Down)] # Put data back into original time resolution

# Settings
numstim <- 1
len <- length(Run)
t <- seq(1, T, by = TR / Down)
tlen <- length(t)

K <- 8
norder <- 4

# Create design matrix
dspline <- Fit_Spline(tc, TR, Runc, T)
Run <- Runc

library(fda)
basisd <- create.bspline.basis(c(0, tlen), nbasis = K + 3, norder = norder)
Bd <- eval.basis(evalarg = t, basisobj = basisd)
Bd <- Bd[, 3:(ncol(Bd) - 1)]

Wi <- matrix(0, nrow = len, ncol = K)
Wji <- tor_make_deconv_mtx3(Run, tlen, 1)$DX
Wi[, 1:K] <- Wji[, 1:tlen] %*% Bd

X <- cbind(1, Wi)
Xd <- X[seq(1, nrow(X), by = Down), ]
X<- X[seq(1, nrow(X), by = Down), ]

# Fit model
b <- ginv(Xd) %*% tc
e <- tc - Xd %*% b
fit <- Xd %*% b
b2 <- b[2:(K + 1)]
hrfd <- Bd %*% b[2:(K + 1)]

# Get parameters
param <- matrix(0, nrow = 3, ncol = numstim)

for (i in 1:numstim) {
  tmp <- get_parameters2(hrfd[, i], 1:length(t))
  param[, i] <- tmp
}

plot(seq(1, 30, by = TR / Down), hrfd)
