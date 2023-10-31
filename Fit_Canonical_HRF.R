Fit_Canonical_HRF <- function(tc, TR, Run, T, p) {
  library(MASS)
  # Fits GLM using canonical hrf (with option of using time and dispersion derivatives)
  
  # INPUTS:
  # tc    - time course
  # TR    - time resolution
  # Runs  - experimental design
  # T     - length of estimated HRF
  # p     - Model type
  # Options: p=1 - only canonical HRF
  #          p=2 - canonical + temporal derivative
  #          p=3 - canonical + time and dispersion derivative
  
  # OUTPUTS:()
  # hrf   - estimated hemodynamic response function
  # fit   - estimated time course
  # e     - residual time course
  # param - estimated amplitude, height and width
  # info  - list containing design matrices, beta values, etc.
  
  d <- length(Run)
  len <- length(Run[[1]])
  t <- seq(1, T, by = TR)
  
  X <- matrix(0, nrow = len, ncol = p * d)
  param <- matrix(0, nrow = 3, ncol = d)
  
  hrf <- list()
  fit <- matrix(0, nrow = len, ncol = d)
  e <- matrix(0, nrow = len, ncol = d)
  
  for (i in 1:d) {
    h <- conv(Run[[i]], CanonicalBasisSet(TR)$h)
    
    X[, (i - 1) * p + 1] <- h[1:len]
    
    if (p > 1) {
      dh <- conv(Run[[i]], CanonicalBasisSet(TR)$dh)
      X[, (i - 1) * p + 2] <- dh[1:len]
    }
    
    if (p > 2) {
      dh2 <- conv(Run[[i]], CanonicalBasisSet(TR)$dh2)
      X[, (i - 1) * p + 3] <- dh2[1:len]
    }
  }
  
  X <- cbind(matrix(1, nrow = len, ncol = 1), X)
  b <- ginv(X) %*% tc
  e <- tc - X %*% b
  fit <- X %*% b
  
  b <- matrix(b[-1], nrow = d, ncol = p, byrow = TRUE)
  bc <- rep(0, d)
  H <- matrix(0, nrow = len, ncol = p)
  
  for (i in 1:d) {
    if (p == 1) {
      bc[i] <- b[i, 1]
      H <- CanonicalBasisSet(TR)$h
    } else if (p == 2) {
      bc[i] <- sign(b[i, 1]) * sqrt(b[i, 1]^2 + b[i, 2]^2)
      H <- cbind(CanonicalBasisSet(TR)$h, CanonicalBasisSet(TR)$dh)
    } else if (p > 2) {
      bc[i] <- sign(b[i, 1]) * sqrt(b[i, 1]^2 + b[i, 2]^2 + b[i, 3]^2)
      H <- cbind(CanonicalBasisSet(TR)$h, CanonicalBasisSet(TR)$dh, CanonicalBasisSet(TR)$dh2)
    }
  }
  
  for (i in 1:d) {
    hrf[[i]] <- H * b[i, i]
    param[, i] <- get_parameters2(hrf[[i]], 1:length(t))
  }
  
  info <- list()
  info$b <- b
  info$bc <- bc
  info$X <- X
  info$H <- H
  
  return(list(hrf = hrf, fit = fit, e = e, param = param, info = info))
}

CanonicalBasisSet <- function(TR) {
  len <- round(40 / TR)
  
  # Define basis function structure with required parameters
  xBF <- list()
  xBF$dt <- TR
  xBF$length <- len
  xBF$name <- 'hrf (with time and dispersion derivatives)'
  
  # Get basis functions using spm_get_bf function (implemented separately)
  xBF <- spm_get_bf(xBF)
  print(xBF$bf)
  
  # Extract basis functions v1, v2, and v3
  v1 <- xBF$bf[1:len, 1]
  v2 <- xBF$bf[1:len, 2]
  v3 <- xBF$bf[1:len, 3]
  
  h <- v1
  dh <- v2 - (crossprod(v2, v1)[1,1] / sum(v1^2)) * v1
  dh2 <- v3 - (crossprod(v3, v1)[1,1] / sum(v1^2)) * v1 - (crossprod(v3, dh)[1,1] / sum(dh^2)) * dh

  h <- max_normalize(h)
  dh <- max_normalize(dh)
  dh2 <- max_normalize(dh2)
  
  return(list(h = h, dh = dh, dh2 = dh2))
}

max_normalize <- function(x){
  if (max(abs(x)) == 0) {
    return (x)
  } else{
    return(x/max(abs(x)))
  }
}
