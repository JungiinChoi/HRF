# Libraries

library(R.matlab)
library(ggplot2)

# Load time course
data <- readMat("/Users/user/Documents/MATLAB/HRF_Est_Toolbox3/tc_0726.mat")
tc <- data$tc
R <- data$R
Run <- data$Run
xsecs <- data$xsecs
hrf <- data$hrf
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

R = c(3, 21, 56, 65, 109, 126, 163, 171, 216, 232, 269, 
      282, 323, 341, 376, 385, 429, 446, 483, 491, 536, 
      552, 589, 602)
Run = rep(0,640)
Run[R] <- 1
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

plot(1:65, hrf, type= "l")
lines(1:59, h1)

pv = ResidScan(e1, FWHM)$p
#PowLoss1 = PowerLoss(e1, fit1, (len-7) , tc, TR, Runc, alpha)

HRF_df$fit1 <- Logit2_fitted$fit
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

#cat("Power Loss:\n")
#cat(PowLoss1, "\n")

# Fit HRF using FIR-model

# Choose mode (FIR/sFIR)

mode = 1;   # 0 - FIR 
            # 1 - smooth FIR

sFIR_fitted <- Fit_sFIR(tc, TR, Run, T, mode)
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

cat("Summary: sFIR\n")

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

#cat("Power Loss:\n")
#cat(PowLoss1, "\n")


# Fit HRF using FIR-model

# Choose mode (FIR/sFIR)

mode = 1;   # 0 - FIR 
# 1 - smooth FIR

sFIR_fitted <- Fit_sFIR(tc, TR, Run, T, mode)
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

cat("Summary: sFIR\n")

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

#cat("Power Loss:\n")
#cat(PowLoss1, "\n")


# Fit HRF using Canonical HRF + 2 derivatives

p = 1

Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Run, 30, p)
h3 = Canonical_HRF_fitted$hrf
fit3 = Canonical_HRF_fitted$fit
e3 = Canonical_HRF_fitted$e
param = Canonical_HRF_fitted$param

pv = ResidScan(e2, FWHM)$p
