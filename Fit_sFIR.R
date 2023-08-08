Fit_sFIR <- function(tc, TR, Run, T, mode) {
  numstim <- length(Run)
  len <- length(Run[[1]])
  t <- seq(1, T, by = TR)
  tlen <- length(t)
  
  Runs <- matrix(0, nrow = len, ncol = numstim)
  for (i in 1:numstim) {
    Runs[, i] <- Run[[i]]
  }
  
  DX <- tor_make_deconv_mtx3(Runs, tlen, 1)$DX
  
  if (mode == 1) {
    C <- outer(1:tlen, 1:tlen, FUN = "-")
    h <- sqrt(1 / (7 / TR))  # 7 seconds smoothing - ref. Goutte
    v <- 0.1
    sig <- 1
    R <- v * exp(-h / 2 * (C^2))
    RI <- solve(R)
    MRI <- matrix(0, nrow = numstim * tlen + 1, ncol = numstim * tlen + 1)
    for (i in 1:numstim) {
      MRI[((i - 1) * tlen + 1):(i * tlen), ((i - 1) * tlen + 1):(i * tlen)] <- RI
    }
    b <- solve(t(DX) %*% DX + sig^2 * MRI) %*% t(DX) %*% tc
    fit <- DX %*% b
    e <- tc - DX %*% b
  } else if (mode == 0) {
    b <- solve(t(DX) %*% DX) %*% t(DX) %*% tc
    fit <- DX %*% b
    e <- tc - DX %*% b
  }
  
  hrf <- matrix(0, nrow = tlen, ncol = numstim)
  param <- matrix(0, nrow = 3, ncol = numstim)
  
  for (i in 1:numstim) {
    hrf[, i] <- t(b[((i - 1) * tlen + 1):(i * tlen)])
    param[, i] <- get_parameters2(hrf[, i], seq(1, tlen))
  }
  
  return(list(hrf = hrf, fit = fit, e = e, param = param))
}
