# HCP Motor Spatial Modeling

## libraries
library(ggplot2)
library(tidyverse)
library(stats)
library(pracma)
library(mgcv)
library(scatterplot3d)
library(Matrix)
library(rgl)


## Load Data
rm(list = ls())
setwd("/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR")
load("Resfit_list_valid.RData")
load("dat_list_valid.RData")
load("Runc_list.RData")
load("valid_i.RData")

## Use valid data only
V <- sum(valid_i)
Tlen <- dim(dat_list_valid[[1]])[2]
Runc <- Runc_list[[1]]
numstim <- length(Runc_list[[1]])
coords <- t(read.table("/Users/user/Downloads/BrainCoordinates.txt"))[valid_i,]


## MGCV modeling of betas for 10 subjects

### DD: Canonical HRF + Derivative
kbeta <- 2

l <- 10

Resfit1 <- Resfit_list_valid[[l]]
dat1 <- dat_list_valid[[l]]

fitted_beta <- matrix(0, nrow = V, ncol = numstim * kbeta)
for (i in 1:V){
  fitted_beta[i,] <- c(t(Resfit1[[i]]$info$b))
}

Blist <- list(length = (numstim * kbeta))
for (i in 1:(numstim * kbeta)){
  gam_df <- data.frame(beta = fitted_beta[,i], x = coords[,1], 
                       y = coords[,2], z = coords[,3])
  gamfit <- gam(beta ~ s(x,y,z, k = 20), data = gam_df, method = "REML")
  Blist[[i]] <- as.matrix(model.matrix(gamfit))
}


#X <- Resfit1[[1]]$info$X[,-1]
#X_ext <- bdiag(replicate(V, X, simplify = FALSE))

#i_perm <- 1:(numstim * kbeta * V)
#j_perm <- c()
#for (i in 1:V){
#  j_perm <- c(j_perm, (0:11)*V+i)
#}
#perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 

#save(X_ext, perm, file = "mgcv_2.RData")

B <- bdiag(Blist)
PB <- perm %*% B
y_ext <- c(t(dat1))

XtX <- t(X_ext) %*% X_ext
BPXtX <- t(PB) %*% XtX
A <- BPXtX %*% PB
Xty_ext <- t(X_ext) %*% y_ext
b <- t(PB) %*% Xty_ext

#gamma_fitted <- matrix(0, ncol = 240, nrow = 10)
gamma_fitted[l,] <- Matrix::solve(A, b)[,1]

#beta_spatial <- array(0, dim = c(10, V, numstim * kbeta))
for (i in 1:(numstim * kbeta)){
  beta_spatial[l,,i] <- (Blist[[i]] %*% gamma_fitted[l,(1+(i-1)*20):(i*20)])[,1]
}

# Nonspatial Beta
Nonspatial_beta <- array(0, dim = c(10, V, numstim * kbeta))
for (l in 1:10){
  for (i in 1:V){
    Nonspatial_beta[l,i,] <- c(t(Resfit_list_valid[[l]][[i]]$info$b))
  }
}

# Fitted Amplitude, Time-to-peak, and Width

A_DD_fit <- T_DD_fit <- W_DD_fit <- array(0, dim = c(10, V, numstim))
for (j in 1:10){
  for (i in 1:V){
    A_DD_fit[j,i, ] <- Resfit_list_valid[[j]][[i]]$param[1,]
    T_DD_fit[j,i, ] <- Resfit_list_valid[[j]][[i]]$param[2,]
    W_DD_fit[j,i, ] <- Resfit_list_valid[[j]][[i]]$param[3,]
  }
}


# Save fitted gamma,beta, and ATW
gamma_fitted_1 <- gamma_fitted
beta_spatial_1 <- beta_spatial
Nonspatial_beta_1 <- Nonspatial_beta
A_DD_fit_1 <- A_DD_fit
T_DD_fit_1 <- T_DD_fit
W_DD_fit_1 <- W_DD_fit
save(gamma_fitted_1, beta_spatial_1, Nonspatial_beta_1, 
     A_DD_fit_1, T_DD_fit_1, W_DD_fit_1, file = "DD_1.RData")


###############################################################
load("DD_1.RData")
### Plot average of spatial betas among 10 subjects

beta_mean <- apply(beta_spatial_1, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
  for (i in (2*j-1):(2*j)){
    cols <- myColorRamp(c("white", "blue"), beta_mean[,i])
    scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
  }
}

### Plot average of nonspatial betas among 10 subjects

Nonspatial_beta <- array(0, dim = c(3, V, numstim * kbeta))
for (l in 1:3){
  for (i in 1:V){
    Nonspatial_beta[l,i,] <- c(t(Resfit_list_valid[[l]][[i]]$info$b))
  }
}

nonspat_beta_mean <- apply(Nonspatial_beta, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(6,3,5)){
  for (i in (2*j-1):(2*j)){
    cols <- myColorRamp(c("white", "blue"), nonspat_beta_mean[,i])
    scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE)
    #legend.col(col = rev(sort(cols)), lev = nonspat_beta_mean[,i])
  }
  title(stim_name[j], line = -3, outer = TRUE)
}

###############################################################

# Plot Non-Spatial Amplitude

A_DD_mean <- apply(A_DD_fit, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
    cols <- myColorRamp(c("white", "blue"), A_DD_mean[,j])
    scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
    #legend.col(col = rev(sort(cols)), lev = beta_mean[,i])
}

# Plot Spatial Amplitude

A_DD_spatial <- array(0, dim = c(10, V, numstim))

TR <- 0.72
H <- cbind(CanonicalBasisSet(TR)$h, CanonicalBasisSet(TR)$dh)

for (i in 1:10){
  for (k in 1:V){
    b <- matrix(beta_spatial[i,k,], nrow = numstim, ncol = 2, byrow = TRUE)
    for (j in 1:numstim) {
      hrf_tmp <- H %*% t(matrix(b[j, ], nrow = 1))
      A_DD_spatial[i,k,j] <- get_parameters2(hrf_tmp, 1:41)[1]
    }
  }
}

A_DD_spat_mean <- apply(A_DD_spatial, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
  cols <- myColorRamp(c("white", "blue"), A_DD_spat_mean[,j])
  scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
  #legend.col(col = rev(sort(cols)), lev = beta_mean[,i])
}

par(mfrow=c(1,1))
boxplot(A_DD_spat_mean, main = "Average Fitted Amplitude")


# Inference on Amplitude

pvalue_A_spat <- matrix(0,nrow = V, ncol = 6)

for (j in 1:6){
  for (i in 1:V){
    pvalue_A_spat[i,j] = t.test(A_DD_spatial[,i,j], alternative = "greater")$p.value
  }
}

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (i in c(6,1,2:5)){
    cols <- ifelse(pvalue_A_spat[,i] < 0.01, "red", "lightgrey")
    scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE)
}


###############################################################

# Inference on Betas

pvalue_beta_spat <- matrix(0,nrow = V, ncol = 12)
pvalue_beta_nspat <- matrix(0,nrow = V, ncol = 12)

for (j in 1:12){
  for (i in 1:V){
    pvalue_beta_spat[i,j] = t.test(beta_spatial[,i,j], alternative = "greater")$p.value
    pvalue_beta_nspat[i,j] = t.test(Nonspatial_beta[,i,j], alternative = "greater")$p.value
  }
}

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(6,3,5)){
  for (i in (2*j-1):(2*j)){
    cols <- ifelse(pvalue_beta_spat[,i] < 0.01, "red", "lightgrey")
    #cols <- myColorRamp(c("red", "white"), pvalue_beta_spat[,i], pvalue_beta_spat[,i])
    scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE)
  }
  #legend.col(col = rev(sort(cols)), lev = pvalue_beta_spat[,i-1])
  title(stim_name[j], line = -3, outer = TRUE)
}


stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(6,3,5)){
  for (i in (2*j-1):(2*j)){
    cols <- ifelse(pvalue_beta_nspat[,i] < 0.11, "red", "lightgrey")
    scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE)
  }
  #legend.col(col = rev(sort(cols)), lev = pvalue_beta_nspat[,i-1])
  title(stim_name[j], line = -3, outer = TRUE)
}


## Fitted non-spatial errors

hrf_nspat <- array(0, dim = c(10, V, Tlen))
e_nspat <- array(0, dim = c(10, V, Tlen))
param_nspat <- array(0, dim = c(10, V, 18))
for(l in 1:10){
  for (i in 1:V){
    Resfitli <- Resfit_list_valid[[l]][[i]]
    hrf_nspat[l,i,] <- Resfitli$fit
    e_nspat[l,i,] <- Resfitli$e
    param_nspat[l,i,] <- Resfitli$param
  }
}

e_nspat_mean <- apply(e_nspat, 3, mean)
e_nspat_sd <- apply(e_nspat, 3, sd)

data = data.frame(t = 1:Tlen, y = e_nspat_mean, sd = e_nspat_sd)
ggplot(data = data, aes(x = t)) + 
  geom_line(aes(y = y), size = 1) + 
  geom_ribbon(aes(y = y, ymin = y-sd, ymax = y+sd), alpha = .2) +
  xlab("Time") + 
  ylab("Timecourse Error")

## Fitted spatial errors

e_spat <- array(0, dim = c(10, V, Tlen))

for(l in 1:10){
  for (i in 1:V){
    betali <- beta_spatial[l,i,]
    e_spat[l,i,] <- (dat_list_valid[[l]][i,] - X %*% betali)[,1]
  }
}

e_spat_mean <- apply(e_spat, 3, mean)
e_spat_sd <- apply(e_spat, 3, sd)

data = data.frame(t = 1:Tlen, y = e_spat_mean, sd = e_spat_sd)
ggplot(data = data, aes(x = t)) + 
  geom_line(aes(y = y), size = 1) + 
  geom_ribbon(aes(y = y, ymin = y-sd, ymax = y+sd), alpha = .2) +
  xlab("Time") + 
  ylab("Timecourse Error")

### in one plot

df_both <- data.frame(t = rep(1:Tlen,2), mean = c(e_nspat_mean, e_spat_mean),
                      sd = c(e_nspat_sd, e_spat_sd), 
                      method = rep(c("non-spatial", "spatial"),each = Tlen))
ggplot(data = df_both, aes(x = t, group = method)) + 
  geom_line(aes(y = mean, color = method), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = method), alpha = .2) +
  xlab("t") +
  ylab("Error") +
  theme_bw() +  
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

boxplot(e_spat_mean, main = "Non-spatial")
boxplot(e_nspat_mean, main = "Spatial")

# Mis-modeling

residscan_nspat <- matrix(0,10,V)
residscan_spat <- matrix(0,10,V)
for (i in 1:10){
  for (j in 1:V){
    residscan_nspat[i,j] <- ResidScan(e_spat[i,j,], 4)$p
    residscan_spat[i,j] <- ResidScan(e_nspat[i,j,], 4)$p
  }
}

residscan_nspat_mean <- apply(residscan_nspat, 2, mean)
residscan_spat_mean <- apply(residscan_spat, 2, mean)

par(mfrow=c(1,2))
cols <- ifelse(residscan_nspat_mean > 0.1, "lightgray", "red")
#cols <- myColorRamp(c("red", "white"), residscan_nspat_mean, c(0,1))
scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
              xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = "Non-spatial")

cols <- ifelse(residscan_spat_mean > 0.1, "lightgray", "red")
#cols <- myColorRamp(c("red", "white"), residscan_spat_mean, c(0,1))

scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
              xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = "Spatial")

legend("right", legend=c("p-values > 0.1", "p-values < 0.1"),
       fill=c("lightgrey", "red"))

boxplot(residscan_nspat_mean, main = "Non-spatial")
boxplot(residscan_spat_mean, main = "Spatial")

##################################

#Spline Model
load("Resfit_list_valid_spline.RData")


## MGCV modeling of betas for 10 subjects

### DD: Canonical HRF + Derivative
kbeta <- 3

gamma_fitted_spline <- matrix(0, ncol = 360, nrow = 10)
beta_spatial_spline <- array(0, dim = c(10, V, numstim * kbeta))
#load("mgcv_2_spline.RData")

for (l in 6:10){
  Resfit1 <- Resfit_list_valid_spline[[l]]
  dat1 <- dat_list_valid[[l]]
  
  fitted_beta_spline <- matrix(0, nrow = V, ncol = numstim * kbeta)
  for (i in 1:V){
    fitted_beta_spline[i,] <- c(t(Resfit1[[i]]$b2))
  }
  Blist <- list(length = (numstim * kbeta))
  for (i in 1:(numstim * kbeta)){
    gam_df <- data.frame(beta = fitted_beta_spline[,i], x = coords[,1], 
                         y = coords[,2], z = coords[,3])
    gamfit <- gam(beta ~ s(x,y,z, k = 20), data = gam_df, method = "REML")
    Blist[[i]] <- as.matrix(model.matrix(gamfit))
  }
  B <- bdiag(Blist)
  PB <- perm %*% B
  y_ext <- c(t(dat1))
  
  XtX <- t(X_ext) %*% X_ext
  BPXtX <- t(PB) %*% XtX
  A <- BPXtX %*% PB
  Xty_ext <- t(X_ext) %*% y_ext
  b <- t(PB) %*% Xty_ext
  
  gamma_fitted_spline[l,] <- Matrix::solve(A, b)[,1]
  
  for (i in 1:(numstim * kbeta)){
    beta_spatial_spline[l,,i] <- (Blist[[i]] %*% gamma_fitted_spline[l,(1+(i-1)*20):(i*20)])[,1]
  }
}

#X <- fit$X[,-1]
#X_ext <- bdiag(replicate(V, X, simplify = FALSE))

#i_perm <- 1:(numstim * kbeta * V)
#j_perm <- c()
#for (i in 1:V){
#  j_perm <- c(j_perm, (0:17)*V+i)
#}
#perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 

#save(X_ext, perm, file = "mgcv_2_spline.RData")

#############################################

hrf_spline <- array(0,dim=c(10,8888,6,41))

for (k in 1:10){
  for (i in 1:8888){
    for (j in 1:6){
      # Get parameters
      hrf_tmp <- fit$B %*% matrix(beta_spatial_spline[k,i,], nrow = 3, ncol = 6)
      hrf_spline[k,i,j,] <- hrf_tmp[,j]
    }
  }
}
hrf_spline_mean <- apply(hrf_spline, c(1,3,4), mean)
hrf_spline_long <- c()
for (i in 1:10){
  for (j in 1:6){
    hrf_spline_long <- c(hrf_spline_long, hrf_spline_mean[i,j,])
  }
}

hrf_spline_df <- data.frame(t = rep(1:41, 10*6), hrf = hrf_spline_long, 
                            stimuli = as.factor(rep(rep(1:6, each = 41),10)),
                            ppl = as.factor(rep(1:10, each = 6*41)))

ggplot(hrf_spline_df, aes(t, hrf, group = interaction(stimuli,ppl))) +
  geom_line(aes(color = stimuli)) +
  labs(title = "Non-spatial HRF Distribution",
       x = "Time",
       y = "y",
       fill = "Stimuli",
       color = "Stimuli") +
  theme_bw() +
  theme(legend.position = "right")


# Non-spatial Amplitude

A_spline_fit <- T_spline_fit <- W_spline_fit <- array(0, dim = c(10, V, numstim))
for (j in 1:10){
  for (i in 1:V){
    A_spline_fit[j,i, ] <- Resfit_list_valid_spline[[j]][[i]]$param[1,]
    T_spline_fit[j,i, ] <- Resfit_list_valid_spline[[j]][[i]]$param[2,]
    W_spline_fit[j,i, ] <- Resfit_list_valid_spline[[j]][[i]]$param[3,]
  }
}
# Plot Non-Spatial Amplitude

A_spline_mean <- apply(A_spline_fit, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
  cols <- myColorRamp(c("white", "blue"), A_spline_mean[,j])
  scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

T_spline_mean <- apply(T_spline_fit, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
  cols <- myColorRamp(c("white", "red"), T_spline_mean[,j])
  scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

W_spline_mean <- apply(W_spline_fit, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
  cols <- myColorRamp(c("white", "green"), W_spline_mean[,j])
  scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

# Plot Spatial Amplitude

A_spline_spatial <- T_spline_spatial <- W_spline_spatial <- array(0, dim = c(10, V, numstim))

B <- fit$B

for (i in 1:10){
  for (k in 1:V){
    b2 <- matrix(beta_spatial_spline[i,k,], nrow = numstim, ncol = 3, byrow = TRUE)
    for (j in 1:numstim) {
      hrf_tmp <- B %*% t(matrix(b2[j,], nrow=1))
      param_tmp <- get_parameters2(hrf_tmp, 1:41)
      A_spline_spatial[i,k,j] <- param_tmp[1]
    }
  }
}

W_spline_spat_mean <- apply(W_spline_spatial, c(2,3), mean)

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
  cols <- myColorRamp(c("white", "green"), W_spline_spat_mean[,j])
  scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

par(mfrow=c(1,1))
cols <- myColorRamp(c("red", "blue"), A_spline_spat_mean[,2] - A_spline_spat_mean[,3])
scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
              xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = "LH - RH")

#####################################################

# Inference on Spatial Betas

p_value_stim <- array(0, dim = c(5, V, numstim))

## F test for nested models.

A_spline_fit <- T_spline_fit <- W_spline_fit <- array(0, dim = c(10, V, numstim))
for (j in 1:10){
  for (i in 1:V){
    A_spline_fit[j,i, ] <- Resfit_list_valid_spline[[j]][[i]]$param[1,]
    T_spline_fit[j,i, ] <- Resfit_list_valid_spline[[j]][[i]]$param[2,]
    W_spline_fit[j,i, ] <- Resfit_list_valid_spline[[j]][[i]]$param[3,]
  }
}

####################################################

# Mis-modeling

e_nspat_spline <- array(0, dim = c(10, V, Tlen))
for(l in 1:10){
  for (i in 1:V){
    Resfitli <- Resfit_list_valid_spline[[l]][[i]]
    e_nspat_spline[l,i,] <- Resfitli$e
  }
}

e_spat <- array(0, dim = c(5, V, Tlen))

for(l in 1:5){
  for (i in 1:V){
    betali <- beta_spatial_spline[l,i,]
    e_spat[l,i,] <- ((dat_list_valid[[l]][i,] - mean(dat_list_valid[[l]][i,])) - X_spline %*% betali)[,1]
  }
}

residscan_nspat <- matrix(0,10,V)
residscan_spat <- matrix(0,5,V)
for (i in 1:5){
  for (j in 1:V){
    residscan_nspat[i,j] <- ResidScan(e_nspat[i,j,], 4)$p
    residscan_spat[i,j] <- ResidScan(e_spat[i,j,], 4)$p
  }
}

for (i in 6:10){
  for (j in 1:V){
    residscan_nspat[i,j] <- ResidScan(e_nspat[i,j,], 4)$p
  }
}

residscan_nspat_mean <- apply(residscan_nspat, 2, mean)
residscan_spat_mean <- apply(residscan_spat, 2, mean)

par(mfrow=c(1,2))
cols <- ifelse(residscan_nspat_mean > 0.1, "lightgray", "red")
scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
              xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = "Non-spatial")

cols <- ifelse(residscan_spat_mean > 0.1, "lightgray", "red")
scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
              xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = "Spatial")

legend("right", legend=c("p-values > 0.1", "p-values < 0.1"),
       fill=c("lightgrey", "red"))

boxplot(residscan_nspat_mean, main = "Non-spatial")
boxplot(residscan_spat_mean, main = "Spatial")

##################################################

#sFIR Model
load("Resfit_list_valid_sFIR.RData")

Resfit_valid_sFIR <- Resfit_list_valid_sFIR[[1]]

# Non-spatial Amplitude

A_sFIR_fit <- T_sFIR_fit <- W_sFIR_fit <- array(0, dim = c(V, numstim))

for (i in 1:V){
  A_sFIR_fit[i, ] <- Resfit_valid_sFIR[[i]]$param[1,]
  T_sFIR_fit[i, ] <- Resfit_valid_sFIR[[i]]$param[2,]
  W_sFIR_fit[i, ] <- Resfit_valid_sFIR[[i]]$param[3,]
}

# Plot Non-Spatial Amplitude


stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")
par(mfrow=c(3,2))
for (j in c(1,6,2:5)){
  cols <- myColorRamp(c("white", "green"), W_sFIR_fit[,j])
  scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}
par(mfrow=c(1,1))
boxplot(W_sFIR_fit)
T_spline_mean <- apply(T_spline_fit, c(2,3), mean)

## MGCV modeling of betas for 10 subjects

### DD: Canonical HRF + Derivative
kbeta <- 41

gamma_fitted_spline <- matrix(0, ncol = 4920, nrow = 10)
beta_spatial_spline <- array(0, dim = c(V, numstim * kbeta))

Resfit1 <- Resfit_list_valid_sFIR[[1]]
dat1 <- dat_list_valid[[1]]
  
fitted_beta_sfir <- matrix(0, nrow = V, ncol = numstim * kbeta)
for (i in 1:V){
  fitted_beta_sfir[i,] <- c(t(Resfit1[[i]]$b[-247,]))
}
  
for (i in 1:1){
  gam_df <- data.frame(beta = fitted_beta_sfir[,i], x = coords[,1], 
                       y = coords[,2], z = coords[,3])
  gamfit <- gam(beta ~ s(x,y,z, k = 20), data = gam_df, method = "REML")
  Blist <- as.matrix(model.matrix(gamfit))
}

B <- bdiag(replicate(246, Blist,  simplify = FALSE))
PB <- perm %*% B
y_ext <- c(t(dat1))
  
XPB <- X_ext %*% PB
XtX <- t(X_ext) %*% X_ext
BPXtX <- crossprod(B,PXtX) 

PXtX <- t(perm) %*% XtX
PXtXP <- PXtX %*% perm
A <- t(B) %*% PXtXP %*% B
A <- BPXtX %*% PB
Xty_ext <- t(X_ext) %*% y_ext
b <- t(PB) %*% Xty_ext
  
gamma_fitted_spline[l,] <- Matrix::solve(A, b)[,1]

for (i in 1:(numstim * kbeta)){
  beta_spatial_spline[l,,i] <- (Blist[[i]] %*% gamma_fitted_spline[l,(1+(i-1)*20):(i*20)])[,1]
}

save(PB, y_ext, BPXtX, file = "sfir.RData")


#X <- Resfit_valid_sFIR[[1]]$DX[,-247]
#X_ext <- bdiag(replicate(V, X, simplify = FALSE))

#i_perm <- 1:(numstim * kbeta*V)
#j_perm <- c()
#for (i in 1:V){
#  j_perm <- c(j_perm, (0:245)*V+i)
#}
#perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 

save(X_ext, perm, file = "mgcv_2_sfir.RData")

