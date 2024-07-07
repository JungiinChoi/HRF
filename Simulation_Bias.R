# Downsampling simulation

# libraries
library(ggplot2)
library(tidyverse)
library(stats)
library(fda)
library(Splinets)

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

# Settings
TR = 1
T = 40
t = seq(1, T, by = TR) 

h <- 1
true_sig <- conv(Run, true_HRF)[1:640]
plot(1:40, true_HRF, type= "l")
plot(1:640, true_sig, type= "l")
tc_noise <- noise_arp(n = 640, phi = c(0.3, 0))
tc <- true_sig + 0.5 * tc_noise
tc <- (tc - mean(tc)) / sd(tc)
xsecs <- 0:40
Run <- rep(0, 640)
Run[100] <- 1
len <- length(Run)

tlen <- length(tc)

# Create design matrix

basis <- create.bspline.basis(rangeval = c(0, tlen), nbasis = 4, norder = 4)
B <- eval.basis(evalarg = 1:tlen, basisobj = basis)

# QR decomposition
qr_result <- qr(B)

# Extract Q and R matrices
Q <- qr.Q(qr_result)
true_beta <- c(0.5,-1,3,-1)
true_HRF <- Q %*% true_beta

plot(1:640, true_HRF, type = "l")

Wi <- matrix(0, nrow = len, ncol = 2)
Wji <- tor_make_deconv_mtx3(matrix(Run, ncol = 1), tlen, 1)$DX
Wi[, 1:2] <- Wji[, 1:tlen] %*% Q

X <- cbind(1, Wi)

# Fit model
b <- ginv(X) %*% tc
b2 <- matrix(b[-1], nrow = 2, ncol = 1)

# Get parameters
hrf <- Q %*% b2
param <- matrix(0, nrow = 3, ncol = numstim)

plot(1:640, hrf, type= "l")

