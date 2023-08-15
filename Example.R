# Libraries

library(R.matlab)
library(ggplot2)

# Load time course
#data <- readMat("/Users/user/Documents/MATLAB/HRF_Est_Toolbox3/tc_0726.mat")
#tc <- data$tc
#R <- data$R
#Run <- data$Run
#xsecs <- data$xsecs
#hrf <- data$hrf
#len <- length(tc)

data <- readMat("/Users/user/Documents/MATLAB/HRF_Est_Toolbox3/tc.mat")
tc <- data$tc
tc <- scale(tc)[,1]
R <- data$R
Run <- data$Run
len <- length(tc)

# Plot of true HRF
ggplot(data = data.frame(xsecs = xsecs[1,], hrf = hrf[,1]), aes(x=xsecs, y=hrf)) +
  geom_line() + 
  theme_minimal() + 
  labs(title="True hrf", x ="t (sec)", y = "")


# Settings

TR = 0.5
T = 30
t = seq(1, T, by = TR)    # samples at which to get Logit HRF Estimate
FWHM = 4                  # FWHM for residual scan
pval = 0.01
df = 600
alpha = 0.001

# Create stick function (sample event onsets)
# Variable R contains onset times
# Variable Run contains stick (a.k.a. delta or indicator) function

#R = c(3, 21, 56, 65, 109, 126, 163, 171, 216, 232, 269, 
      282, 323, 341, 376, 385, 429, 446, 483, 491, 536, 
      552, 589, 602)
#Run = rep(0,640)
#Run[R] <- 1
Runc <- list(Run)
onset <- rep(-100,640)
onset[R] <- R

HRF_df <- data.frame(t = 1:len, tc = tc, onset = onset)
  
ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc, color = 'Data')) +
  theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 640)) +
  scale_color_manual(name = "", values=c('blue'))


# Fit HRF using IL-function
# Choose mode (deterministic/stochastic)

mode = 0    # 0 - deterministic approach 
            # 1 - simulated annealing approach
            # Please note that when using simulated annealing approach you
            # may need to perform some tuning before use.

Logit2_fitted <- Fit_Logit2(tc, TR, Runc, T, mode)
h1 = Logit2_fitted$hrf
fit1 = Logit2_fitted$fit
e1 = Logit2_fitted$e
param = Logit2_fitted$param

pv = ResidScan(e1, FWHM)$p

HRF_df$fit1 <- fit1
ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc, color = 'Data')) +
  geom_line(aes(x=t, y=fit1, color = 'IL')) +
  theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 640)) +
  scale_color_manual(name = "", values=c('blue', 'red'))

cat("Summary: IL_function\n")

cat("Amplitude:\n")
cat(param[1], "\n")

cat("Time-to-peak:\n")
cat(param[2] * TR, "\n")

cat("Width:\n")
cat(param[3] * TR, "\n")

cat("MSE:\n")
cat((1 / (len - 1) * sum(e1^2)), "\n")

cat("Mis-modeling:\n")
cat(pv, "\n")

# Fit HRF using FIR-model

# Choose mode (FIR/sFIR)

mode = 1;   # 0 - FIR 
            # 1 - smooth FIR

sFIR_fitted <- Fit_sFIR(tc, TR, Runc, T, mode)
h2 = sFIR_fitted$hrf
fit2 = sFIR_fitted$fit
e2 = sFIR_fitted$e
param = sFIR_fitted$param

pv = ResidScan(e2, FWHM)$p

#[PowLoss2] = PowerLoss(e2, fit2, (len-T) , tc, TR, Runc, alpha);

HRF_df$fit2 <- fit2
ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc, color = 'Data')) +
  geom_line(aes(x=t, y=fit1, color = 'IL')) +
  geom_line(aes(x=t, y=fit2, color = 'sFIR')) +
  theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 640)) +
  scale_color_manual(name = "", values=c('blue', 'red', 'green'))

plot(1:65, hrf, type = "line", xlab = "", ylim = c(-0.5,1.15))
lines(1:59, h1, , col = "red")
lines(1:59, h2, type = "line", col = "green")

cat("Summary: sFIR\n")

cat("Amplitude:\n")
cat(param[1], "\n")

cat("Time-to-peak:\n")
cat(param[2] * TR, "\n")

cat("Width:\n")
cat(param[3] * TR, "\n")

cat("MSE:\n")
cat((1 / (len - 1) * sum(e2^2)), "\n")

cat("Mis-modeling:\n")
cat(pv, "\n")


# Fit HRF using Canonical HRF + 2 derivatives

p = 1

Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Runc, 30, p)
h3 = Canonical_HRF_fitted$hrf[[1]]
fit3 = Canonical_HRF_fitted$fit[, 1]
e3 = Canonical_HRF_fitted$e[, 1]
param = Canonical_HRF_fitted$param[,1]

pv = ResidScan(e2, FWHM)$p

HRF_df$fit3 <- fit3
ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc, color = 'Data')) +
  geom_line(aes(x=t, y=fit1, color = 'IL')) +
  geom_line(aes(x=t, y=fit2, color = 'sFIR')) +
  geom_line(aes(x=t, y=fit3, color = 'DD')) +
  theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 640)) +
  scale_color_manual(name = "", values=c('blue', 'magenta', 'red', 'green'))

plot(1:65, hrf, type = "line", xlab = "", ylim = c(-0.5,1.15))
lines(1:59, h1, col = "red")
lines(1:59, h2, col = "green")
lines(1:60, h3, col = "magenta")

cat('Summary: Canonical + 2 derivatives\n')
cat('Amplitude\n', param[1], '\n')
cat('Time-to-peak\n', param[2] * TR, '\n')
cat('Width\n', param[3] * TR, '\n')

cat('MSE:\n', (1 / (len - 1) * sum(e3^2)), '\n')
cat('Mis-modeling\n', pv, '\n')

# Fit HRF using B-splies

Spline_fitted <- Fit_Spline(tc, TR, Runc, 30)
h4 <- Spline_fitted$hrf
fit4 <- Spline_fitted$fit
e4 <- Spline_fitted$e
param <- Spline_fitted$param
pv <- ResidScan(e4, FWHM)$p

HRF_df$fit4 <- fit4
ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc, color = 'Data')) +
  geom_line(aes(x=t, y=fit1, color = 'IL')) +
  geom_line(aes(x=t, y=fit2, color = 'sFIR')) +
  geom_line(aes(x=t, y=fit3, color = 'DD')) +
  geom_line(aes(x=t, y=fit4, color = 'Spline')) +
  theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 640)) +
  scale_color_manual(name = "", values=c('blue', 'magenta', 'red', 'green', 'skyblue'))


cat("Summary: B-spline\n")
cat("Amplitude: ", param[1], "\n")
cat("Time-to-peak: ", param[2] * TR, "\n")
cat("Width: ", param[3] * TR, "\n")
cat("MSE: ", (1 / (len - 1) * sum(e4^2)), "\n")
cat("Mis-modeling: ", pv, "\n")



# Fit HRF using non-linear gamma function

NLgamma_fitted <- Fit_NLgamma(tc, TR, Runc, 30)
h5 = NLgamma_fitted$hrf[,1]
fit5 = NLgamma_fitted$fit
e5 = NLgamma_fitted$e
param = NLgamma_fitted$param[,1]
pv = ResidScan(e5, FWHM)$p

HRF_df$fit5 <- fit5
ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc, color = 'Data')) +
  geom_line(aes(x=t, y=fit1, color = 'IL')) +
  geom_line(aes(x=t, y=fit2, color = 'sFIR')) +
  geom_line(aes(x=t, y=fit3, color = 'DD')) +
  geom_line(aes(x=t, y=fit4, color = 'Spline')) +
  geom_line(aes(x=t, y=fit5, color = 'NL')) +
  theme_minimal() +
  labs(title = "Sample time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -4, xend = onset, yend = -3)) +
  coord_cartesian(xlim = c(0, 640)) +
  scale_color_manual(name = "", values=c('blue', 'magenta', 'red', 'yellow', 'green', 'skyblue'))

cat('Summary: Non-linear gamma\n')
cat('Amplitude\n', param[1], '\n')
cat('Time-to-peak\n', param[2] * TR, '\n')
cat('Width\n', param[3] * TR, '\n')
cat('MSE:\n', (1 / (len - 1) * sum(e5^2)), '\n')
cat('Mis-modeling\n', pv, '\n')


# HRF 
plot(1:65, rep(0,65), type = "line", xlab = "", ylim = c(-0.5,1.15))
lines(1:59, h1, col = "red")
lines(1:59, h2, col = "green")
lines(1:59, h3[-1], col = "magenta")
lines(1:59, h4, col = "skyblue")
lines(1:59, h5, col = "yellow")



