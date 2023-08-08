library(stats)

spm_hrf <- function(RT, p = NULL, T = NULL) {
  # Default parameters
  P <- c(6, 16, 1, 1, 6, 0, 32)
  
  # Use default parameters if P is not provided
  if (!is.null(p)) {
    P[1:length(p)] <- p
  }
  
  # Microtime resolution
  if (is.null(T)) {
    T <- 16
  }
  
  # Modelled haemodynamic response function - {mixture of Gammas}
  dt <- RT / T
  u <- 0:(ceiling(P[7] / dt)) - P[6] / dt
  hrf <- (dgamma(u, shape = P[1] / P[3], scale = dt / P[3]) -
            dgamma(u, shape = P[2] / P[4], scale = dt / P[4]) / P[5])
  hrf <- hrf[(1:(floor(P[7] / RT))) * T]
  hrf <- hrf / sum(hrf)
  
  return(list(hrf = hrf, p = P))
}



# Probability Density Function (PDF) of Gamma distribution
spm_Gpdf <- function(x, h, l){
  return(dgamma(x, shape = h, scale = l))
}

# Function to get basis functions for event-related responses in R
spm_get_bf <- function(xBF) {
  dt <- xBF$dt
  
  # Model event-related responses
  if (!"name" %in% names(xBF)) {
    message("Hemodynamic Basis functions...")
    Ctype <- c(
      'hrf',
      'hrf (with time derivative)',
      'hrf (with time and dispersion derivatives)',
      'Fourier set',
      'Fourier set (Hanning)',
      'Gamma functions',
      'Finite Impulse Response'
    )
    str <- "Select basis set"
    Sel <- utils::menu(Ctype)
    xBF$name <- Ctype[Sel]
  }
  
  # Get order and length parameters
  if (xBF$name %in% c('Fourier set', 'Fourier set (Hanning)', 
                      'Finite Impulse Response', 'Gamma functions')){
    l <- if (exists("xBF$length")) xBF$length else utils::readline("window length {secs}", 32)
    xBF$length <- l
    
    h <- if (exists("xBF$order")) xBF$order else utils::readline("order", 4)
    xBF$order <- h
  } 
  
  # Create basis functions
  if (xBF$name == 'Fourier set' || xBF$name == 'Fourier set (Hanning)') {
    pst <- seq(0, l, dt)
    pst <- pst / max(pst)
    g <- if (xBF$name == 'Fourier set (Hanning)') (1 - cos(2 * pi * pst)) / 2 else rep(1, length(pst))
    bf <- g
    for (i in 1:h) {
      bf <- cbind(bf, g * sin(i * 2 * pi * pst))
      bf <- cbind(bf, g * cos(i * 2 * pi * pst))
    }
  } else if (xBF$name == 'Gamma functions') {
    pst <- seq(0, l, dt)
    bf <- sapply(2:(1 + h), function(i) {
      m <- 2^i
      s <- sqrt(m)
      spm_Gpdf(pst, (m/s)^2, m/s^2)
    })
  } else if (xBF$name == 'Finite Impulse Response') {
    bin <- l / h
    bf <- kronecker(diag(h), matrix(1, round(bin/dt), 1))
  } else if (xBF$name == 'NONE') {
    bf <- 1
  } else {
    fMRI_T <- if (exists("xBF$T")) xBF$T else 16
    bfp <- spm_hrf(dt, NULL, fMRI_T)
    bf <- bfp$hrf
    bf <- matrix(bf, ncol = 1)
    p <- bfp$p
    print("first bf")
    print(bf)
    
    if (grepl("time", xBF$name)) {
      dp <- 1
      
      bf <- cbind(bf, (bf[, 1] - spm_hrf(dt, c(p[1:5], p[6] + dp), fMRI_T)$hrf) / dp)
      
      if (grepl("dispersion", xBF$name)) {
        dp <- 0.01
        bf <- cbind(bf, (bf[, 1] - spm_hrf(dt, c(p[1:2], p[3] + dp, p[4:6]), fMRI_T)$hrf) / dp)
      }
    }
    xBF$length <- nrow(bf) * dt
    xBF$order <- ncol(bf)
  }
  
  # Orthogonalise and fill in basis function structure
  print("bf")
  print(bf)
  
  xBF$bf <- qr(bf)$qr
  
  return(xBF)
}

spm_orth <- function(X, OPT = "pad") {
  # Recursive Gram-Schmidt orthogonalization
  n <- nrow(X)
  m <- ncol(X)
  X <- X[, colSums(X != 0) > 0]  # Remove zero columns
  
  rankX <- qr(X)$rank
  x <- X[, 1]
  j <- 1
  
  for (i in 2:m) {
    D <- X[, i]
    D <- D - x %*% (solve(crossprod(x)) %*% D)
    if (norm(D, type = "1") > exp(-32)) {
      x <- cbind(x, D)
      j <- c(j, i)
    }
    if (length(j) == rankX) break
  }
  
  # Normalization, if requested
  if (OPT == "pad") {
    X <- matrix(0, n, m)
    X[, j] <- x
  } else if (OPT == "norm") {
    X <- x / sqrt(sum(x^2))
  } else {
    X <- x
  }
  return(X)
}
