library(fda)

K <- 8                  # Number of b-spline basis
norder <- 4             # Order of b-spline basis

Fit_Spline_Sim <- function(tc, TR, Run, T, K = 8, norder = 4) {
  numstim <- length(Run)   # Number conditions
  len <- length(Run[[1]])  # length of run
  t <- seq(from = 1, to = T, by = TR)
  tlen <- length(t)        # Number of time points in HRF
  
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

h <- 1

true_sig <- b * conv(Run, hrfmat[h,])[1:640]

fitted_hrf_spline <- array(0, dim = c(7,4,40,100))
params_spline <- array(0, dim = c(7,4,3,100))
MSE_spline <- array(0, dim = c(7,4,2,100))

for (i in 1:100){
  tc_noise <- noise_arp(n = 640, phi = c(0.3, 0))
  tc <- true_sig + 0.5 * tc_noise
  tc <- (tc - mean(tc)) / sd(tc)
  xsecs <- 0:40
  
  for (K in 5:11){
    for (norder in 3:6){
      Spline_fitted <- Fit_Spline_Sim(tc, TR, Runc, T, K, norder)
      fitted_hrf_spline[K-4, norder-2, ,i] <- Spline_fitted$hrf
      params_spline[K-4, norder-2, ,i] <- Spline_fitted$param
      e4 <- Spline_fitted$e
    }
  }
}

bias_spline_A <- matrix(0, nrow = 7, ncol = 4)
bias_spline_T <- matrix(0, nrow = 7, ncol = 4)
bias_spline_W <- matrix(0, nrow = 7, ncol = 4)

for (K in 5:11){
  for (norder in 3:6){
    bias_spline_A[K-4, norder-2] <- mean(params_spline[K-4, norder-2, 1,] - true_params[h,1])
    bias_spline_T[K-4, norder-2] <- mean(params_spline[K-4, norder-2, 2,] - true_params[h,2])
    bias_spline_W[K-4, norder-2] <- mean(params_spline[K-4, norder-2, 3,] - true_params[h,3])
  }
}

Spline_W <- data.frame(K = rep(5:11,each = 4), norder = rep(3:6, 7), bias = abs(c(t(bias_spline_W))))

ggplot(Spline_W, aes(x = factor(K), y = factor(norder), fill = bias)) +
  geom_tile(color = "black") +
  coord_fixed() +
  xlab("Number of b-spline basis")+
  ylab("Order of b-spline basis") +
  scale_fill_gradient(low = "white", high = "green", limits = c(0, 2)) 



abs(scale(bias_spline_A) + scale(bias_spline_T) + scale(bias_spline_W)) 


