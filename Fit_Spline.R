library(fda)

Fit_Spline <- function(tc, TR, Run, T) {
  numstim <- length(Run)   # Number conditions
  len <- length(Run[[1]])  # length of run
  t <- seq(from = 1, to = T, by = TR)
  tlen <- length(t)        # Number of time points in HRF
  
  K <- 8                  # Number of b-spline basis
  norder <- 4             # Order of b-spline basis
  
  # Create design matrix
  basis <- create.bspline.basis(rangeval = c(0, tlen), nbasis = K + 3, norder = norder)
  B <- eval.basis(evalarg = 1:tlen, basisobj = basis)
  B <- B[, 3:(ncol(B) - 1)]
  
  Wi <- matrix(0, nrow = len, ncol = numstim * K)
  for (j in 1:numstim) {
    Wji <- tor_make_deconv_mtx3(Run[[j]], tlen, 1)$DX
    Wi[, ((j - 1) * K + 1):(j * K)] <- Wji[, 1:tlen] %*% B
    
  }
  
  X <- cbind(1, Wi)
  
  # Fit model
  b <- ginv(X) %*% tc
  e <- tc - X %*% b
  fit <- X %*% b
  
  b2 <- matrix(b[-1], nrow = K, ncol = numstim)
  
  # Get parameters
  hrf <- B %*% b2
  param <- matrix(0, nrow = 3, ncol = numstim)
  
  for (i in 1:numstim) {
    tmp <- get_parameters2(hrf[, i], 1:length(t))
    print("tmp")
    print(tmp)
    print(param)
    param[, i] <- tmp
  }
  
  return(list(hrf = hrf, fit = fit, e = e, param = param))
}
