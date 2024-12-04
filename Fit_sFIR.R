Fit_sFIR <- function(tc, TR, Run, T, mode) {
  # Fits FIR and smooth FIR model
  #
  # INPUTS:
  #
  # tc    - time course
  # TR    - time resolution
  # Run   - experimental design
  # T     - length of estimated HRF
  # mode  - FIR or smooth FIR
  #   options:
  #       0 - standard FIR 
  #       1 - smooth FIR
  #
  # OUTPUTS:
  #
  # hrf   - estimated hemodynamic response function
  # fit   - estimated time course
  # e     - residual time course
  # param - estimated amplitude, height and width
  #
  # Created by Martin Lindquist on 10/02/09
  # Last edited: 05/17/13 (ML)
  
  numstim <- length(Run)
  len <- length(Run[[1]])
  t <- seq(1, T, TR)
  tlen <- length(t)
  
  Runs <- matrix(0, len, numstim)
  for (i in 1:numstim) {
    Runs[, i] <- Run[[i]]
  }
  
  DX <- tor_make_deconv_mtx3(Runs, tlen, 1)$DX
  
  if (mode == 1) {
    C <- matrix(rep(t, each = tlen), nrow = tlen)
    h <- sqrt(1 / (7 / TR))  # 7 seconds smoothing - ref. Goutte
    
    v <- 0.1
    sig <- 1
    
    R <- v * exp(-h / 2 * (C - t(C))^2)
    RI <- base::solve(R)
    MRI <- matrix(0, numstim * tlen + 1, numstim * tlen + 1)
    for (i in 1:numstim) {
      MRI[((i - 1) * tlen + 1):(i * tlen), ((i - 1) * tlen + 1):(i * tlen)] <- RI
    }
    b <- base::solve(t(DX) %*% DX + sig^2 * MRI) %*% t(DX) %*% tc
    fit <- DX %*% b
    e <- tc - DX %*% b
    
  } else if (mode == 0) {
    b <- MASS::ginv(DX) %*% tc
    fit <- DX %*% b
    e <- tc - DX %*% b
  }
  
  hrf <- matrix(0, tlen, numstim)
  param <- matrix(0, 3, numstim)
  
  for (i in 1:numstim) {
    hrf[, i] <- b[((i - 1) * tlen + 1):(i * tlen)]
    param[, i] <- get_parameters2(hrf[, i], seq(1, tlen))
  }
  
  list(hrf = hrf, fit = fit, e = e, param = param, DX = DX, b = b)
}
