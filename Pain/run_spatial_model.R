# Load Rscripts

source("../Fit_Canonical_HRF.R")
source("../Fit_Spline.R")
source("../get_parameters2.R")
source("../spm_hrf.R")
source("../tor_make_deconv_mtx3.R")

library(mgcv)
library(Matrix)

# Run Nonspatial Model

load("data_list.RData")
load("Runc.RData")
load("reglen.RData")
numreg <- 268
TR <- 0.460
numstim <- 4

Nonspatial_beta_DD <- array(0, dim = c(numreg, max(reglen), numstim * 2))
Nonspatial_beta_spline <- array(0, dim = c(numreg, max(reglen), numstim * 8))
A_DD_fit <- T_DD_fit <- W_DD_fit <- A_spline_fit <- T_spline_fit <- W_spline_fit <- array(0, dim = c(numreg,max(reglen), numstim))

e_DD_fit <- e_sp_fit <- array(0, dim = c(numreg, max(reglen), 872))
                              
for (k in 1:numreg){
  dat <- dat_list[[k]] 
  coord <- coord_list[[k]]
  V <- nrow(dat)
  numstim <- 4
  
  # Nonspatial Beta
  for (i in 1:V){
    fit_DD <- Fit_Canonical_HRF(dat[i,], TR, Runc, 30, 2)
    fit_spline <- Fit_Spline(dat[i,], TR, Runc, 30)
    e_DD_fit[k,i,] <- fit_DD$e
    e_sp_fit[k,i,] <- fit_spline$e
  }
  print(k)
}

save(Nonspatial_beta_DD, Nonspatial_beta_spline, A_DD_fit, A_spline_fit, 
     file = "/Users/user/Documents/JHU/research/HRF_Est/Pain/nonspatial.RData")

save(e_DD_fit, e_sp_fit,
     file = "/Users/user/Documents/JHU/research/HRF_Est/Pain/residuals.RData")

# Spatial Thin Plate Regression Model

## MGCV modeling of betas for 57 regions

### DD: Canonical HRF + Derivative

beta_spatial <- list()
gamma_fitted <- matrix(0, nrow = numreg, ncol = 8 * 20)

for (k in 1:numreg){
  coords <- coord_list[[k]]
  dat1 <- dat_list[[k]] 
  V <- nrow(dat1)

  gam_df <- data.frame(beta = Nonspatial_beta_DD[k,1:V,1], x = coords[,1], 
                       y = coords[,2], z = coords[,3])
  gamfit <- gam(beta ~ s(x,y,z, k = 20), data = gam_df, method = "REML")
  Blist <- as.matrix(model.matrix(gamfit))

  X <- fit_DD$info$X[,-1]
  X_ext <- bdiag(replicate(V, X, simplify = FALSE))
  
  i_perm <- 1:(numstim * 2 * V)
  j_perm <- c()
  for (i in 1:V){
    j_perm <- c(j_perm, (0:(2*4-1))*V+i)
  }
  perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 
  
  B <- bdiag(replicate(numstim * 2, Blist, simplify = FALSE))
  PB <- perm %*% B
  y_ext <- c(t(dat1))
  
  XtX <- t(X_ext) %*% X_ext
  BPXtX <- t(PB) %*% XtX
  A <- BPXtX %*% PB
  Xty_ext <- t(X_ext) %*% y_ext
  b <- t(PB) %*% Xty_ext
  
  gamma_fitted[k,] <- Matrix::solve(A, b)[,1]
  beta_spatial[[k]] <- matrix(0, V, numstim * 2)
  for (i in 1:(numstim * 2)){
    beta_spatial[[k]][,i] <- (Blist %*% gamma_fitted[k,(1+(i-1)*20):(i*20)])[,1]
  }
  print(k)
}


### Spatial Amplitude

A_DD_spatial <- array(0, dim = c(numreg, max(reglen), numstim))

H <- cbind(CanonicalBasisSet(TR)$h, CanonicalBasisSet(TR)$dh)

for (n in 1:numreg){
  for (i in 1:reglen[n]){
    b <- matrix(beta_spatial[[n]][i,], nrow = numstim, ncol = 2, byrow = TRUE)
    for (j in 1:numstim) {
      hrf_tmp <- H %*% t(matrix(b[j, ], nrow = 1))
      A_DD_spatial[n,i,j] <- get_parameters2(hrf_tmp, 1:length(seq(1, 30, by = TR)))[1]
    }
  }
}


### Spline: B-spline


beta_spatial_sp <- list()
gamma_fitted_sp <- matrix(0, nrow = numreg, ncol = 20 * 4 * 8)

X <- fit_spline$X[,-1]

for (k in 1:numreg){
  coords <- coord_list[[k]]
  dat1 <- dat_list[[k]] 
  V <- nrow(dat1)
  
  gam_df <- data.frame(beta = Nonspatial_beta_spline[k,1:V,1], x = coords[,1], 
                       y = coords[,2], z = coords[,3])
  gamfit <- gam(beta ~ s(x,y,z, k = 20), data = gam_df, method = "REML")
  Blist <- as.matrix(model.matrix(gamfit))
  
  X_ext <- bdiag(replicate(V, X, simplify = FALSE))
  XtX <- t(X_ext) %*% X_ext
  
  i_perm <- 1:(numstim * 8 * V)
  j_perm <- c()
  for (i in 1:V){
    j_perm <- c(j_perm, (0:(4 * 8 - 1))*V+i)
  }
  perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 

  B <- bdiag(replicate(numstim * 8, Blist, simplify = FALSE))
  PB <- perm %*% B
  y_ext <- c(t(dat1))
  
  BPXtX <- t(PB) %*% XtX
  A <- BPXtX %*% PB
  Xty_ext <- t(X_ext) %*% y_ext
  b <- t(PB) %*% Xty_ext
  
  gamma_fitted_sp[k,] <- Matrix::solve(A, b)[,1]
  beta_spatial_sp[[k]] <- matrix(0, V, numstim * 8)
  for (i in 1:(numstim * 8)){
    beta_spatial_sp[[k]][,i] <- (Blist %*% gamma_fitted_sp[k,(1+(i-1)*20):(i*20)])[,1]
  }
  print(k)
}

### Spatial Amplitude

A_sp_spatial <- array(0, dim = c(numreg, max(reglen), numstim))

B <- fit_spline$B

for (n in 1:numreg){
  for (i in 1:reglen[n]){
    b2 <- matrix(beta_spatial_sp[[n]][i,], nrow = numstim, ncol = 8, byrow = TRUE)
    for (j in 1:numstim) {
      hrf_tmp <- B %*% t(matrix(b2[j,], nrow=1))
      A_sp_spatial[n,i,j] <- get_parameters2(hrf_tmp, 1:length(seq(1, 30, by = TR)))[1]
    }
  }
}

