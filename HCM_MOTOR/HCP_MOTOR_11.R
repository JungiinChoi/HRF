# The first number in each triplet is the onset (in seconds) of the period
# the second number is the duration (in seconds) of the period, and 
# the third number is the value of the input during that period. 

# libraries
library(ggplot2)
library(tidyverse)
library(stats)
library(pracma)

TR <- 0.72
nvol <- 284

folders <- list.dirs("/Users/user/Personal Research Dropbox/Choi Jungin/MOTOR_DATA", 
                          recursive = FALSE)
K <- 10

dat_list <- vector("list", length = K)
Runc_list <- vector("list", length = K)

# Save timecourse data and stimuli function for each person

for (j in 1:1){
  setwd(folders[j])
  
  # LR 
  X1_LR <- matrix(0, ncol = 3, nrow = nvol)
  X1_LR[,1] <- 1
  X1_LR[,2] <- (1:nvol)/nvol
  X1_LR[,3] <- ((1:nvol)^2)/(nvol^2)
  M_LR <- read.table('tfMRI_MOTOR_LR/Movement_Regressors.txt')
  X_LR <- as.matrix(cbind(X1_LR, M_LR))
  
  # fmri timecourse data for LR
  dat1 <- read.table('tfMRI_MOTOR_LR/tfMRI_MOTOR_LR_SomMotNetwork.txt', sep = ",",
                     header = FALSE, quote = "", stringsAsFactors = FALSE)
  
  # timecourse error term (- projection on regressors)
  for (i in 1:nrow(dat1)) {
    y <- t(dat1[i,])
    e <- y - X_LR %*% solve(t(X_LR) %*% X_LR, t(X_LR)) %*% y
    dat1[i,] <- t(e)
  }
  
  # cue
  cue <- round(read.table('tfMRI_MOTOR_LR/EVs/cue.txt')/TR)
  S1_LR <- rep(0, nvol)
  for (i in 1:nrow(cue)) {
    S1_LR[cue[i,1]:(cue[i,1] + cue[i,2])] <- 1
  }
  
  # left foot
  lf <- round(read.table('tfMRI_MOTOR_LR/EVs/lf.txt')/TR)
  S2_LR <- rep(0, nvol)
  for (i in 1:nrow(lf)) {
    S2_LR[lf[i,1]:(lf[i,1] + lf[i,2])] <- 1
  }
  
  # right foot
  rf <- round(read.table('tfMRI_MOTOR_LR/EVs/rf.txt')/TR)
  S3_LR <- rep(0, nvol)
  for (i in 1:nrow(rf)) {
    S3_LR[rf[i,1]:(rf[i,1] + rf[i,2])] <- 1
  }
  
  # left hand
  lh <- round(read.table('tfMRI_MOTOR_LR/EVs/lh.txt')/TR)
  S4_LR <- rep(0, nvol)
  for (i in 1:nrow(lh)) {
    S4_LR[lh[i,1]:(lh[i,1] + lh[i,2])] <- 1
  }
  
  # right hand
  rh <- round(read.table('tfMRI_MOTOR_LR/EVs/rh.txt')/TR)
  S5_LR <- rep(0, nvol)
  for (i in 1:nrow(rh)) {
    S5_LR[rh[i,1]:(rh[i,1] + rh[i,2])] <- 1
  }
  
  # tongue
  ton <- round(read.table('tfMRI_MOTOR_LR/EVs/t.txt')/TR)
  S6_LR <- rep(0, nvol)
  for (i in 1:nrow(ton)) {
    S6_LR[ton[i,1]:(ton[i,1] + ton[i,2])] <- 1
  }
  
  # RL
  X2_RL <- matrix(0, ncol = 3, nrow = nvol)
  X2_RL[,1] <- 1
  X2_RL[,2] <- (1:nvol)/nvol
  X2_RL[,3] <- ((1:nvol)^2)/(nvol^2)
  M_RL <- read.table('tfMRI_MOTOR_RL/Movement_Regressors.txt')
  X_RL <- as.matrix(cbind(X2_RL, M_RL))
  
  # fmri timecourse data for RL
  dat2 <- read.table('tfMRI_MOTOR_RL/tfMRI_MOTOR_RL_SomMotNetwork.txt', sep = ",",
                     header = FALSE, quote = "", stringsAsFactors = FALSE)
  
  # timecourse error term (- projection on regressors)
  for (i in 1:nrow(dat2)) {
    y <- t(dat2[i,])
    e <- y - X_RL %*% solve(t(X_RL) %*% X_RL, t(X_RL)) %*% y
    dat2[i,] <- t(e)
  }
  
  # cue
  cue <- round(read.table('tfMRI_MOTOR_RL/EVs/cue.txt')/TR)
  S1_RL <- rep(0, nvol)
  for (i in 1:nrow(cue)) {
    S1_RL[cue[i,1]:(cue[i,1] + cue[i,2])] <- 1
  }
  
  # left foot
  lf <- round(read.table('tfMRI_MOTOR_RL/EVs/lf.txt')/TR)
  S2_RL <- rep(0, nvol)
  for (i in 1:nrow(lf)) {
    S2_RL[lf[i,1]:(lf[i,1] + lf[i,2])] <- 1
  }
  
  # right foot
  rf <- round(read.table('tfMRI_MOTOR_RL/EVs/rf.txt')/TR)
  S3_RL <- rep(0, nvol)
  for (i in 1:nrow(rf)) {
    S3_RL[rf[i,1]:(rf[i,1] + rf[i,2])] <- 1
  }
  
  # left hand
  lh <- round(read.table('tfMRI_MOTOR_RL/EVs/lh.txt')/TR)
  S4_RL <- rep(0, nvol)
  for (i in 1:nrow(lh)) {
    S4_RL[lh[i,1]:(lh[i,1] + lh[i,2])] <- 1
  }
  
  # right hand
  rh <- round(read.table('tfMRI_MOTOR_RL/EVs/rh.txt')/TR)
  S5_RL <- rep(0, nvol)
  for (i in 1:nrow(rh)) {
    S5_RL[rh[i,1]:(rh[i,1] + rh[i,2])] <- 1
  }
  
  # tongue
  ton <- round(read.table('tfMRI_MOTOR_RL/EVs/t.txt')/TR)
  S6_RL <- rep(0, nvol)
  for (i in 1:nrow(ton)) {
    S6_RL[ton[i,1]:(ton[i,1] + ton[i,2])] <- 1
  }
  
  # Experimental Designs
  Runc <- list(
    c(S1_LR, S1_RL),
    c(S2_LR, S2_RL),
    c(S3_LR, S3_RL),
    c(S4_LR, S4_RL),
    c(S5_LR, S5_RL),
    c(S6_LR, S6_RL)
  )
  
  # Delete initial 4 data for RL and LR
  dat <- as.matrix(cbind(dat1, dat2))[,-c(1:4,nvol+(1:4))]
  for (i in 1:length(Runc)){
    Runc[[i]] <- Runc[[i]][-c(1:4,nvol+(1:4))]
  }
  
  dat_list[[j-10]] <- dat
  Runc_list[[j-10]] <- Runc
}

for (k in 1:10){
  Runc <- Runc_list[[k]]
  for (i in 1:length(Runc)){
    Runc[[i]] <- Runc[[i]][-c(1:4,nvol+(1:4))]
  }
  Runc_list[[k]] <- Runc
}


num <- 9374
setwd("/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR")
save(dat_list, file = "dat_list.RData")
save(Runc_list, file = "Runc_list.RData")

# Plot Stimuli Function: Same for all patients
Runc <- Runc_list[[1]]
df <- data.frame(
  x = rep(1:length(Runc[[1]]), length(Runc)),
  y = unlist(Runc),
  S = factor(rep(c("cue", "lf", "rf", "lh", "rh", "tongue"), each = length(Runc[[1]])))
)

ggplot(df, aes(x = x, y = y, group = S)) +
  geom_line(aes(color = S)) +
  labs(title = "Stimuli function ", x = "t", y = "") 


# Debugging timecourse data
V <- nrow(dat_list[[1]])
Tlen <- ncol(dat_list[[1]])

sum0 <- rep(0,V)
for (i in 1:K){
  sum0 <- sum0 + apply(dat_list[[i]], 1, function(x) {sum(x == 0) >= Tlen/2})
}

valid_i <- sum0 <=5
setwd("/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR")
save(valid_i, file = "valid_i.RData")

coords <- t(read.table("/Users/user/Downloads/BrainCoordinates.txt"))

scatterplot3d(x = coords[,1], y = coords[,2], z = coords[,3],
              xlab = "x", ylab = "y", zlab = "z", pch = 16, color = ifelse(valid_i, "lightgray", "red"),
              main = "Invalid Voxels")

# Plot average timecourse for 10 subjects
df_tc <- data.frame(t = rep(1:Tlen, 10), subject = as.factor(rep(1:10, each = Tlen)),
                    mean = c(sapply(dat_list, function(x) {apply(x[valid_i,],2,mean)})),
                    sd = c(sapply(dat_list, function(x) {apply(x[valid_i,],2,sd)})))

df_tc_tmp <- df_tc 
ggplot(data = df_tc_tmp, aes(x = t, group = subject)) + 
  geom_line(aes(y = mean, color = subject), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = subject), alpha = .1) +
  xlab("t") + 
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  title("Distribution of Timecourse Data")

## Let's use only valid voxels
# Canonical DD method for 10 ppl

Resfit_list_valid <- vector("list", length = K)
V <- sum(valid_i)
for (k in 11:20){
  Resfit <- list()
  dat <- dat_list[[k-10]][valid_i,]
  Runc <- Runc_list[[1]]
  for (i in 1:V) {
    fit <- Fit_Canonical_HRF(dat[i,], TR, Runc, 30, 2)
    Resfit[[i]] <- fit
  }
  Resfit_list_valid[[k-10]] <- Resfit
}

Runc_list_valid <- Runc_list

dat_list_valid <- lapply(dat_list, function(x){x[valid_i,]})
setwd("/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR")
save(Resfit_list_valid, file = "Resfit_list_valid.RData")
save(dat_list_valid, file = "dat_list_valid.RData")
save(Runc_list_valid, file = "Runc_list_valid.RData")


#########################################

# B-Spline method for 10 ppl

Resfit_list_valid_spline <- vector("list", length = K)
for (k in 1:10){
  Resfit <- list()
  dat <- dat_list_valid[[k]]
  Runc <- Runc_list[[1]]
  for (i in 1:nrow(dat_list_valid[[1]])){
    fit <- Fit_Spline(dat[i,], TR, Runc, 30)
    Resfit[[i]] <- fit
  }
  Resfit_list_valid_spline[[k]] <- Resfit
}

## Plot HRF
plot(1:41, fit$hrf[,1], type = "l", col = "white", ylim = c(-30,100))
Resfit<- Resfit_list_valid_spline[[1]] 
for (i in 1:V){
  lines(1:41, Resfit[[i]]$hrf[,5])
}

## Plot basis of HRF
plot(1:41,B[,3])

setwd("/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR")
save(Resfit_list_valid_spline, file = "Resfit_list_valid_spline.RData")

#########################################

# sFIR method for 10 ppl

Resfit_list_valid_sFIR <- vector("list", length = 10)
V <- nrow(dat_list_valid[[1]])
for (k in 1:1){
  Resfit <- list()
  dat <- dat_list_valid[[k]]
  Runc <- Runc_list[[1]]
  for (i in 1:V) {
    fit <- Fit_sFIR(dat[i,], TR, Runc, 30, 1)
    Resfit[[i]] <- fit
  }
  Resfit_list_valid_sFIR[[k]] <- Resfit
}


for (k in 1:1){
  dat <- dat_list_valid[[k]]
  Runc <- Runc_list[[1]]
  for (i in 1:1) {
    fit <- Fit_sFIR(dat[i,], TR, Runc, 30, 1)
    DX <- fit$DX
    MRI <- fit$MRI
  }
}

par(mfrow=c(1,1))
plot(1:560, DX[,244], type="l")
for (i in 2:247){
  lines(1:560,DX[,i])
}

setwd("/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR")
save(Resfit_list_valid_sFIR, file = "Resfit_list_valid_sFIR.RData")


## Plot Estimated Amplitude, time-to-peak, and Width: Check Spatial Dependencies

coord <- read.table("/Users/user/Downloads/BrainCoordinates.txt")
coord <- t(coord)
library(scatterplot3d)
library(graphics)
library(stats)

## For the first patient

A_DD_fit <- T_DD_fit <- W_DD_fit <- matrix(0, nrow = num, ncol = 6)
for (i in 1:num){
  A_DD_fit[i, ] <- Resfit_list[[1]][[i]]$param[1,]
  T_DD_fit[i, ] <- Resfit_list[[1]][[i]]$param[2,]
  W_DD_fit[i, ] <- Resfit_list[[1]][[i]]$param[3,]
}

beta_fitted <- matrix(0, nrow = num, ncol = 12)
for (i in 1:num){
  beta_fitted[i,] <- c(t(Resfit_list[[1]][[i]]$info$b))
}

library(rgl)

myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

myColorRamp_minmax <- function(colors, values, min, max) {
  v <- (values - min)/(max-min)
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

cols <- myColorRamp(c("white", "black"), beta_fitted[,1]) 

### Amplitude
plot3d(x = coord[,1], y = coord[,2], z = coord[,3], col = cols,
       xlab = xlab, ylab = ylab, zlab = zlab, size = 7)
scatterplot3d(coord[,1],coord[,2],coord[,3], pch = 16, color=cols,
              xlab = xlab, ylab = ylab, zlab = zlab)

## For all patients (average)

A_DD_fit_mean <- T_DD_fit_mean <- W_DD_fit_mean <- matrix(0, nrow = num, ncol = 6)
for (k in 1:10){
  for (i in 1:num){
    A_DD_fit_mean[i, ] <- A_DD_fit_mean[i, ] + Resfit_list[[k]][[i]]$param[1,]
    T_DD_fit_mean[i, ] <- T_DD_fit_mean[i, ] + Resfit_list[[k]][[i]]$param[2,]
    W_DD_fit_mean[i, ] <- W_DD_fit_mean[i, ] + Resfit_list[[k]][[i]]$param[3,]
  }
}


### For the first patient

beta <- array(0, dim=c(num, 6, 2))
p_value <- array(0, dim=c(num, 6, 2))
p_value_stim <- matrix(0, nrow = num, ncol = 6)

Resfit <- Resfit_list[[1]]

for (i in 1:num){
  Restmp <- Resfit[[i]]
  beta_tmp <- Restmp$info$b
  X <- Restmp$info$X[,-1]
  beta[i,,] <- beta_tmp
  
  lm_df <- data.frame(X)
  lm_df$Y <- dat[i,]
  # Fit a linear regression model
  lm_model <- lm(Y ~ ., data = lm_df)
  
  p_value[i,,] <- matrix(coef(summary(lm_model))[-1, "Pr(>|t|)"], ncol = 2, byrow=T)
  
  beta_hat <- lm_model$coefficients[-1]
  df <- length(beta_hat) - 1
  se_beta_hat <- coef(summary(lm_model))[-1, "Std. Error"]
  
  for (j in 1:6){
    vector <- rep(0,12)
    vector[(2*j-1):(2*j)] <- 1
    
    # Calculate the t-statistic
    t_statistic <- sum(beta_hat * vector) / sqrt(sum((vector * se_beta_hat)^2))
    
    # Calculate the p-value
    p_value_stim[i,j] <- 2 * pt(abs(t_statistic), df = df, lower.tail = FALSE)
  }
}

d <- as.data.frame(p_value_stim)
d[is.na(d)] <- 1

dev.off()
par(mfrow = c(2,3))

for (i in 1:6){
  cols <- myColorRamp(c("darkred", "white"), log(d)[,i])
  scatterplot3d(coord[,1],coord[,2],coord[,3], pch = 16, color=cols,
                xlab = xlab, ylab = ylab, zlab = zlab)
}

d_canonical <- as.data.frame(p_value[,,1])
d_canonical[is.na(d_canonical)] <- 1

dev.off()
par(mfrow = c(2,3))
for (i in 1:6){
  cols <- myColorRamp(c("darkred", "white"), log(d_canonical)[,i])
  scatterplot3d(coord[,1],coord[,2],coord[,3], pch = 16, color=cols,
                xlab = xlab, ylab = ylab, zlab = zlab)
}


# B-Spline for 10 ppl

Spline_fit_list <- list(length = 10)

for (k in 1:10){
  Resfit <- list()
  dat <- dat_list[[k]]
  Runc <- Runc_list[[k]]
  for (i in 1:num) {
    Resfit[[i]] <- Fit_Spline(dat[i,], TR, Runc, 30)
  }
  Spline_fit_list[[k]] <- Resfit
  print(paste0(k," is done"))
}

for (k in 1:1){
  dat <- dat_list_valid[[k]]
  Runc <- Runc_list[[k]]
  for (i in 1:1) {
    X_spline <- Fit_Spline(dat[i,], TR, Runc, 30)$X[,-1]
  }
}
par(mfrow=c(1,1))
plot(1:560,X_spline[,8],type="l")


## compute p-values for 6 null hypothesis

beta_Spline <- array(0, dim=c(num, 6, 5))
p_value_Spline <- array(0, dim=c(num, 6, 5))
p_value_stim_Spline <- matrix(0, nrow = num, ncol = 6)

Resfit <- Spline_fit_list[[1]]
dat <- dat_list[[1]]

for (i in 1:num){
  Restmp <- Resfit[[i]]
  beta_tmp <- t(Restmp$b2)
  beta_Spline[i,,] <- beta_tmp
  X <- Wi
  lm_df <- data.frame(X)
  lm_df$Y <- dat[i,]
  
  # Fit a linear regression model
  lm_model <- lm(Y ~ ., data = lm_df)
  
  p_value_Spline[i,,] <- matrix(coef(summary(lm_model))[-1, "Pr(>|t|)"], ncol = 2, byrow=T)
  
  beta_hat <- lm_model$coefficients[-1]
  df <- length(beta_hat) - 1
  se_beta_hat <- coef(summary(lm_model))[-1, "Std. Error"]
  
  for (j in 1:6){
    vector <- rep(0,30)
    vector[(5*j-4):(5*j)] <- 1
    
    # Calculate the t-statistic
    t_statistic <- sum(beta_hat * vector) / sqrt(sum((vector * se_beta_hat)^2))
    
    # Calculate the p-value
    p_value_stim_Spline[i,j] <- 2 * pt(abs(t_statistic), df = df, lower.tail = FALSE)
  }
}


d_Spline <- as.data.frame(p_value_Spline[,,4])
d_Spline[is.na(d_Spline)] <- 1

dev.off()
par(mfrow = c(2,3))
for (i in 1:6){
  cols <- myColorRamp(c("darkred", "white"), log(d_Spline)[,i])
  scatterplot3d(coord[,1],coord[,2],coord[,3], pch = 16, color=cols,
                xlab = xlab, ylab = ylab, zlab = zlab)
}





## Plot estimated HRF function and Timecourse

## Plot Estimated Amplitude, time-to-peak, and Width: Check Spatial Dependencies


## For the first patient

A_Spline_fit <- T_Spline_fit <- W_Spline_fit <- matrix(0, nrow = num, ncol = 6)
for (i in 1:num){
  A_Spline_fit[i, ] <- Spline_fit_list[[1]][[i]]$param[1,]
  T_Spline_fit[i, ] <- Spline_fit_list[[1]][[i]]$param[2,]
  W_Spline_fit[i, ] <- Spline_fit_list[[1]][[i]]$param[3,]
}

### Amplitude
dev.off()
par(mfrow = c(2,3))
for (i in 1:6){
  cols <- myColorRamp(c("white", "darkred"), A_Spline_fit[,i])
  scatterplot3d(coord[,1],coord[,2],coord[,3], pch = 16, color=cols,
                xlab = xlab, ylab = ylab, zlab = zlab)
}
cols <- myColorRamp(c("white", "darkred"), A_Spline_fit[,2])
plot3d(x = coord[,1], y = coord[,2], z = coord[,3], col = cols,
       xlab = xlab, ylab = ylab, zlab = zlab, size = 7)

### Amplitude
par(mfrow = c(2,3))
for (i in 1:6){
  cols <- myColorRamp(c("white", "darkgreen"), W_Spline_fit[,i])
  scatterplot3d(coord[,1],coord[,2],coord[,3], pch = 16, color=cols,
                xlab = xlab, ylab = ylab, zlab = zlab)
}

