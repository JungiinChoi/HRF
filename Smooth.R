# Run HRF estimation using smoothed timecourse data

smooth_dat <- readMat("/Users/user/Documents/JHU/research/HRF_Est/Pain/smo_6.mat")$smo
smo_dat_6 <- matrix(0,nrow(smooth_dat), ncol(smooth_dat))

for (i in 1:nrow(smo_dat_6)) {
  y <- matrix(smooth_dat[i,], ncol = 1)
  e <- y - X %*% base::solve(t(X) %*% X, t(X)) %*% y
  smo_dat_6[i,] <- t(e)
  print(i)
}


## Load Functions

source("../Fit_Canonical_HRF.R")
source("../Fit_Spline.R")
source("../get_parameters2.R")
source("../spm_hrf.R")
source("../tor_make_deconv_mtx3.R")

library(mgcv)
library(Matrix)

## Run Nonspatial Model

load("Runc.RData")
load("reglen.RData")
numreg <- 268
TR <- 0.460
numstim <- 4
V <- nrow(smo_dat_3)

A_DD_fit <- matrix(0, nrow = V, ncol = 4)

numstim <- 4

for (k in 1:V){
  fit_DD <- Fit_Canonical_HRF(smo_dat_6[k,], TR, Runc, 30, 2)
  A_DD_fit[k,] <- fit_DD$param[1,]
  print(k)
}

write.csv(A_DD_fit, file = "/Users/user/Documents/JHU/research/HRF_Est/Pain/A_DD_fit_smo_6_st.csv", 
          row.names = FALSE)


# Plot smoothed timecourse

e_sp_sec1 <- dat_list[[1]][1:300,]
e_sp_sec1 <- smo_dat_3[1:300,]
e_sp_sec1_mean <- apply(e_sp_sec1,2,mean)
sec1 <- data.frame(t = 1:872, y = c(t(e_sp_sec1)))

ggplot(sec1) + 
  geom_line(aes(x=t, y = y, color = "Data"), alpha=0.2) +
  geom_line(data = data.frame(t = 1:872, y = e_sp_sec1_mean),
            aes(x = t, y = y, color = "Mean")) +
  scale_color_manual(name = "", values = c("Data" = "blue", "Mean" = "black"))+
  ggtitle("Timecourse Data after Smoothing for Region 1") +
  ylim(-2000,2000)+
  theme_bw()

