Fit_NLgamma <- function(tc, TR, Run, T) {
  numstim <- length(Run)
  len <- length(Run[[1]])
  t <- seq(1, T, by = TR)
  
  V0 <- matrix(c(1, 6, 1), nrow = 3, ncol = numstim, byrow = TRUE)  # [Height, Delay, Onset]
  
  # Find optimal values
  VM <- matrix(0, nrow = 3, ncol = numstim)
  
  for (g in seq_len(numstim)) {
    VM[, g] <- optim(V0[, g], msq_nl_gamma, Run = Run, TR = TR, T = T, tc = tc)$par
  }
  
  # Use optimal values to fit hemodynamic response functions
  hrf <- matrix(0, nrow = length(t), ncol = numstim)
  fitt <- matrix(0, nrow = len, ncol = numstim)
  param <- matrix(0, nrow = 3, ncol = numstim)
  
  for (g in seq_len(numstim)) {
    hrf[, g] <- NL_gamma(TR, T, VM[, g])  # Calculate HRF estimate (fit, given theta)
    param[, g] <- get_parameters2(hrf[, g], t)
    fits <- conv(Run[[g]], hrf[, g])
    fitt[, g] <- fits[seq_len(len)]
  }
  
  fit <- rowSums(fitt)
  e <- tc - fit
  
  return(list(hrf = hrf, fit = fit, e = e, param = param))
}

####################################################################
# SUBFUNCTIONS
####################################################################

msq_nl_gamma <- function(V, Run, TR, T, tc) {
  numstim <- length(Run)
  len <- length(Run[[1]])
  t <- seq(1, T, by = TR)
  h <- matrix(0, nrow = length(t), ncol = numstim)
  yhatt <- matrix(0, nrow = len, ncol = numstim)
  
  for (k in seq_len(numstim)) {
    print(NL_gamma(TR, T, V))
    print(length(NL_gamma(TR, T, V)))
    print(length(h[, k]))
    h[, k] <- NL_gamma(TR, T, V)  # Get NL gamma model corresponding to parameters V
    yhat <- conv(Run[[k]], h[, k])  # Convolve NL model with stick function
    yhatt[, k] <- yhat[seq_len(len)]
  }
  
  yhat2 <- rowSums(yhatt)  # Sum models together to get overall estimate
  
  return(sum((tc - yhat2)^2))  # Calculate cost function
}

NL_gamma <- function(TR, T, V) {
  # Inverse logit - creates fitted curve from parameter estimates
  # t: vector of time points
  # V: parameters
  
  # 3 logistic functions to be summed together
  height <- V[1]
  delay <- 6
  udelay <- 16
  dispersion <- 1
  udisp <- 1
  rtou <- 6
  onset <- 1
  klength <- T - 1
  
  normhrf <- spm_hrf(TR, c(delay, udelay, dispersion, udisp, rtou, onset, klength))$hrf
  h <- height * normhrf / max(normhrf)
  
  return(h)
}
