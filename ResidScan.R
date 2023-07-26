ResidScan <- function(res, FWHM) {
  res_ns <- res
  res <- res / sd(res)
  len <- length(res)
  
  # Create Gaussian Kernel
  sig <- ceiling(FWHM / (2 * sqrt(2 * log(2))))
  klen <- 3 * sig
  kern <- dnorm((-klen:klen), mean = 0, sd = sig)
  kern <- kern / sqrt(sum(kern^2))
  
  # Convolve
  x <- stats::convolve(res, kern, type = "open")
  sres <- x[(klen + 1):(length(x) - klen)]
  
  x <- stats::convolve(res_ns, kern / sum(kern), type = "open")
  sres_ns <- x[(klen + 1):(length(x) - klen)]
  
  # Find Max value
  a <- max(abs(sres))
  
  # Find p-values using Gaussian Random Field theory
  z <- Euler_p(1, a, len, FWHM)
  z <- 2 * z # Two-sided test
  p <- min(1, z)
  
  return(list(p = p, sres = sres, sres_ns = sres_ns))
}

# Subfunction
Euler_p <- function(myDim, value, N, fwhm) {
  # Constants
  myfactor <- 4 * log(2)
  pi2 <- 2 * pi
  exptsq <- exp(-(value^2) / 2)
  
  # Euler Characteristic Densties
  rho <- numeric(5)
  rho[1] <- 1 - pnorm(value)
  rho[2] <- myfactor^(0.5) * exptsq / pi2
  rho[3] <- myfactor * exptsq * value / (pi2^(1.5))
  rho[4] <- myfactor^(1.5) * exptsq * (value^2 - 1) / (pi2^2)
  rho[5] <- myfactor^2 * exptsq * (value^3 - 3 * value) / (pi2^(5/2))
  
  # Resel Count
  R0 <- 1
  R1 <- N / fwhm
  
  # P-value
  pval <- R0 * rho[1] + R1 * rho[2]
  
  return(pval)
}
