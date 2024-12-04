# Libraries

library(R.matlab)
library(ggplot2)

R <- readMat("/Users/user/Documents/MATLAB/HRF_Est_Toolbox3/stimuli.mat")$R
Run <- matrix(0, nrow = 640, ncol = 1)
Run[R,1] <- 1
len <- 640

hrf <- readMat("/Users/user/Documents/MATLAB/HRF_Est_Toolbox3/true_hrf.mat")$hrf[,1]
xsecs <- seq(0, 40, by = 0.5)

tc_mat <- matrix(0,nrow = len, ncol = 100)
#fitted_hrf_IL <- matrix(0, nrow = 30, ncol = 100)
#params_IL <- matrix(0, nrow = 3, ncol = 100)
#MSE_IL <- matrix(0, nrow = 2, ncol = 100)
fitted_hrf_sFIR <- matrix(0, nrow = 30, ncol = 100)
params_sFIR <- matrix(0, nrow = 3, ncol = 100)
MSE_sFIR <- matrix(0, nrow = 2, ncol = 100)
fitted_hrf_DD <- matrix(0, nrow = 30, ncol = 100)
params_DD <- matrix(0, nrow = 3, ncol = 100)
MSE_DD <- matrix(0, nrow = 2, ncol = 100)
fitted_hrf_spline <- matrix(0, nrow = 30, ncol = 100)
params_spline <- matrix(0, nrow = 3, ncol = 100)
MSE_spline <- matrix(0, nrow = 2, ncol = 100)
#fitted_hrf_NL <- matrix(0, nrow = 30, ncol = 100)
#params_NL <- matrix(0, nrow = 3, ncol = 100)
#MSE_NL <- matrix(0, nrow = 2, ncol = 100)

# Settings

TR = 1
T = 40
t = seq(1, T, by = TR)    # samples at which to get Logit HRF Estimate
FWHM = 4                  # FWHM for residual scan
pval = 0.01
df = 600
alpha = 0.001

# Create stick function (sample event onsets)
# Variable R contains onset times
# Variable Run contains stick (a.k.a. delta or indicator) function

Runc <- list(Run)
onset <- rep(-100,640)
onset[R] <- R

HRF_df <- data.frame(t = 1:len, tc = tc, onset = onset)
set.seed(906)

for (i in 1:100){
  tc_noise <- noise_arp(n = 640, phi = c(0.3, 0))
  tc <- true_sig + 0.5 * tc_noise
  tc <- (tc - mean(tc)) / sd(tc)
  xsecs <- 0:32
  tc_mat[,i] <- tc
  
  #IL
  #Logit2_fitted <- Fit_Logit2(tc, TR, Runc, T, 0)
  #fitted_hrf_IL[,i] <- Logit2_fitted$hrf
  #params_IL[,i] <- Logit2_fitted$param
  #e1 <- Logit2_fitted$e
  #MSE_IL[,i] <- c((1 / (len - 1) * sum(e1^2)), ResidScan(e1, FWHM)$p)
  
  #sFIR
  sFIR_fitted <- Fit_sFIR(tc, TR, Runc, T, 1)
  fitted_hrf_sFIR[,i] <- sFIR_fitted$hrf
  params_sFIR[,i] <- sFIR_fitted$param
  e2 = sFIR_fitted$e
  MSE_sFIR[,i] <- c((1 / (len - 1) * sum(e2^2)), ResidScan(e2, FWHM)$p)
  
  #Canonical HRF + 2 derivatives
  Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Runc, T, 1)
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
  #NLgamma_fitted <- Fit_NLgamma(tc, TR, Runc, T)
  #fitted_hrf_NL[,i] <- NLgamma_fitted$hrf[,1]
  #params_NL[,i] <- NLgamma_fitted$param[,1]
  #e5 = NLgamma_fitted$e
  #MSE_NL[,i] <- c((1 / (len - 1) * sum(e5^2)), ResidScan(e5, FWHM)$p)
}

p <- ggplot(data.frame(xsecs = xsecs[1:30], hrf = hrfmat[1,])) + 
  geom_line(aes(x=xsecs, y = hrf), color = "black", size = 1.5)

for (i in 1:100){
  #hrf_IL <- fitted_hrf_IL[,i]  
  hrf_sFIR <- fitted_hrf_sFIR[,i]
  hrf_DD <- fitted_hrf_DD[,i]
  hrf_Spline <- fitted_hrf_spline[,i]
  #hrf_NL <- fitted_hrf_NL[,i]
  
  #p <- p + geom_line(data = data.frame(xsecs = xsecs[1:30], hrf = hrf_IL),
  #  aes(x = xsecs, y = hrf, colour = "IL"), alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:30], hrf = hrf_sFIR),
                     aes(x = xsecs, y = hrf, color = "sFIR"), alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:30], hrf = hrf_DD),
                     aes(x = xsecs, y = hrf, color = "DD"), alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:30], hrf = hrf_Spline),
                     aes(x = xsecs, y = hrf, color = "Spline"), alpha=0.2)
  #p <- p + geom_line(data = data.frame(xsecs = xsecs[1:30], hrf = hrf_NL),
  #                   aes(x = xsecs, y = hrf, color = "NL"), alpha=0.2)
}

p + scale_color_manual(name = "", values = c("IL" = "red", "sFIR" = "green", "DD" = "magenta", "Spline" = "blue", "NL" = "yellow"))


p <- ggplot()
  
for (i in 1:100){
  p <- p +   geom_line(data = data.frame(t = 1:640, tc = tc_mat[,i]),
                       aes(x = t, y = tc), colour = "blue", alpha=0.2)
}

true_sig <- scale(readMat("/Users/user/Documents/MATLAB/HRF_Est_Toolbox3/true_sig.mat")$true.sig)[,1]
p <- p +
  geom_line(data = data.frame(t = 1:640, tc = true_sig),
            aes(x = t, y = tc), colour = "black")

p + theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 640))

for (i in 1:100){
  hrf_IL <- fitted_hrf_IL[,i]  
  hrf_sFIR <- fitted_hrf_sFIR[,i]
  hrf_DD <- fitted_hrf_DD[,i]
  hrf_Spline <- fitted_hrf_spline[,i]
  hrf_NL <- fitted_hrf_NL[,i]
  
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:59], hrf = hrf_IL),
                     aes(x = xsecs, y = hrf), colour = "red", alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:59], hrf = hrf_sFIR),
                     aes(x = xsecs, y = hrf), colour = "green", alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:60], hrf = hrf_DD),
                     aes(x = xsecs, y = hrf), colour = "magenta", alpha=0.2)
  p <- p + geom_line(data = data.frame(xsecs = xsecs[1:59], hrf = hrf_NL),
                     aes(x = xsecs, y = hrf), colour = "blue", alpha=0.2)
}


params_df <- params_IL
par(mfrow=c(1,3))
boxplot(params_IL[1,], main = "Amplitude", xlab = paste0("mean: ", round(mean(params_IL[1,]),2), "  sd: ", round(sd(params_IL[1,]),2)))
boxplot(params_IL[2,], main = "Time-to-peak", xlab = paste0("mean: ", round(mean(params_IL[2,]),2), "  sd: ", round(sd(params_IL[2,]),2)))
boxplot(params_IL[3,], main = "Width", xlab = paste0("mean: ", round(mean(params_IL[3,]),2), "  sd: ", round(sd(params_IL[3,]),2)))

boxplot(params_sFIR[1,])
boxplot(params_sFIR[2,], main = "Time-to-peak")
boxplot(params_sFIR[3,], main = "Width")

boxplot(params_DD[1,], main = "Amplitude")
boxplot(params_DD[2,], main = "Time-to-peak")
boxplot(params_DD[3,], main = "Width")

boxplot(params_NL[1,], main = "Amplitude")
boxplot(params_NL[2,], main = "Time-to-peak")
boxplot(params_NL[3,], main = "Width")

params_DD


