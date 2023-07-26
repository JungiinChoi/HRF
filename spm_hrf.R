spm_hrf <- function(RT, P, T = 16) {
  # Parameters of the response function
  if (missing(P)) {
    tryCatch({
      P <- c(6, 16, 1, 1, 6, 0, 32)
    }, error = function(e) {
      stop("Parameters not provided and default values not available.")
    })
  }
  
  if (length(P) > 7) {
    stop("Number of parameters should be at most 7.")
  }
  
  # Microtime resolution
  if (missing(T)) {
    tryCatch({
      T <- 16
    }, error = function(e) {
      stop("Microtime resolution not provided and default value not available.")
    })
  }
  
  # Modelled haemodynamic response function - {mixture of Gammas}
  dt <- RT / T
  u <- seq(0, ceiling(P[7] / dt)) - P[6] / dt
  hrf <- spm_Gpdf(u, P[1] / P[3], dt / P[3]) - spm_Gpdf(u, P[2] / P[4], dt / P[4]) / P[5]
  hrf <- hrf[seq(1, floor(P[7] / RT) * T, T)]
  hrf <- hrf / sum(hrf)
  
  return(list(hrf = hrf, p = P))
}


# Probability Density Function (PDF) of Gamma distribution
spm_Gpdf <- function(x, h, l){
  return(dgamma(x, shape = h, scale = l))
}

# Function to get basis functions for event-related responses in R
spm_get_bf <- function(xBF) {
  # Length of time bin
  if (missing(xBF$dt)) {
    str <- "time bin for basis functions {secs}"
    xBF$dt <- utils::readline(prompt = str, default = 1/16)
  }
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
  switch(xBF$name,
         'Fourier set',
         'Fourier set (Hanning)',
         'Gamma functions',
         'Finite Impulse Response', {
           l <- if (exists("xBF$length")) xBF$length else utils::readline("window length {secs}", 32)
           xBF$length <- l
           
           h <- if (exists("xBF$order")) xBF$order else utils::readline("order", 4)
           xBF$order <- h
         })
  
  # Create basis functions
  switch(xBF$name,
         'Fourier set',
         'Fourier set (Hanning)', {
           pst <- seq(0, l, dt)
           pst <- pst / max(pst)
           g <- if (xBF$name == 'Fourier set (Hanning)') (1 - cos(2 * pi * pst)) / 2 else rep(1, length(pst))
           bf <- g
           for (i in 1:h) {
             bf <- cbind(bf, g * sin(i * 2 * pi * pst))
             bf <- cbind(bf, g * cos(i * 2 * pi * pst))
           }
         },
         'Gamma functions', {
           pst <- seq(0, l, dt)
           bf <- sapply(2:(1 + h), function(i) {
             m <- 2^i
             s <- sqrt(m)
             spm_Gpdf(pst, (m/s)^2, m/s^2)
           })
         },
         'Finite Impulse Response', {
           bin <- l / h
           bf <- kronecker(diag(h), matrix(1, round(bin/dt), 1))
         },
         'NONE', {
           bf <- 1
         },
         {
           fMRI_T <- if (exists("xBF$T")) xBF$T else spm_get_defaults("stats.fmri.t")
           bf <- spm_hrf(dt, NULL, fMRI_T)
           
           if (grepl("time", xBF$name)) {
             dp <- 1
             bf <- cbind(bf, (bf[, 1] - spm_hrf(dt, c(p[1:5], p[6] + dp), fMRI_T)) / dp)
             
             if (grepl("dispersion", xBF$name)) {
               dp <- 0.01
               bf <- cbind(bf, (bf[, 1] - spm_hrf(dt, c(p[1:2], p[3] + dp, p[4:6]), fMRI_T)) / dp)
             }
           }
           xBF$length <- nrow(bf) * dt
           xBF$order <- ncol(bf)
         })
  
  # Orthogonalise and fill in basis function structure
  xBF$bf <- spm_orth(bf)
  
  return(xBF)
}
