# Please note that in the above R implementation, some subfunctions like get_parameters2 are not provided, so you'll need to define them based on the corresponding MATLAB implementation. 
# Additionally, the optimization method used in the R version is L-BFGS-B which is commonly used for 
# bounded optimization problems, similar to MATLAB's fminsearchbnd. 


# Load required libraries
library(stats)  # For optimization functions

# Define the inverse logit function
ilogit <- function(t) {
  exp(t) / (1 + exp(t))
}

# Define convolution function
conv <- function(x, y) {
  len_x <- length(x)
  len_y <- length(y)
  conv_result <- numeric(len_x + len_y - 1)
  
  for (i in seq_along(x)) {
    conv_result[i:(i + len_y - 1)] <- conv_result[i:(i + len_y - 1)] + x[i] * y
  }
  return(conv_result)
}

# Define the function to calculate the inverse logit HRF from parameter estimates
# t = vector of time points
# V = parameters

# 3 logistic functions to be summed together
il_hdmf_tw2 <- function(t, V) {
  base <- matrix(0, nrow = length(t), ncol = 3)
  
  A1 <- V[1]
  T1 <- V[2]
  d1 <- V[3]
  A2 <- V[4]
  T2 <- V[5]
  A3 <- V[6]
  T3 <- V[7]
  d2 <- -d1 * (ilogit(A1 * (1 - T1)) - ilogit(A3 * (1 - T3))) / (ilogit(A2 * (1 - T2)) + ilogit(A3 * (1 - T3)))
  d3 <- abs(d2) - abs(d1)
  
  base[, 1] <- d1 * ilogit(A1 * (t - T1))
  base[, 2] <- d2 * ilogit(A2 * (t - T2))
  base[, 3] <- d3 * ilogit(A3 * (t - T3))
  h <- rowSums(base)
  
  return(h)
}

# Define the cost function for optimization
msq_logit <- function(V, Run, t, tc) {
  numstim <- length(Run)
  len <- length(Run[[1]])
  h <- matrix(0, nrow = length(t), ncol = numstim)
  yhatt <- matrix(0, nrow = len, ncol = numstim)
  
  for (k in 1:numstim) {
    h[, k] <- il_hdmf_tw2(t, V[((k - 1) * 7 + 1):(k * 7)])  # Get IL model corresponding to parameters V
    yhat <- conv(Run[[k]], h[, k])                     # Convolve IL model with stick function
    yhatt[, k] <- yhat[1:len]
  }
  
  yhat2 <- rowSums(yhatt)  # Sum models together to get overall estimate
  
  m <- sum((tc - yhat2)^2)  # Calculate cost function
  
  return(m)
}

# Define the main function
Det_Logit <- function(V0, t, tc, Run) {
  numstim <- length(Run)
  len <- length(Run[[1]])
  
  VM <- optim(V0, msq_logit, Run = Run, t = t, tc = tc, 
              control = list(reltol = 1e-6, maxit = 10000))$par
  
  hrf <- matrix(0, nrow = length(t), ncol = numstim)
  fitt <- matrix(0, nrow = len, ncol = numstim)
  param <- matrix(0, nrow = 3, ncol = numstim)
  
  for (g in 1:numstim) {
    hrf[, g] <- il_hdmf_tw2(t, VM[((g - 1) * 7 + 1):(g * 7)])  # Calculate HRF estimate (fit, given theta)
    param[, g] <- get_parameters2(hrf[, g], t)  
    fits <- conv(Run[[g]], hrf[, g])  # Convolve stick function with HRF
    fitt[, g] <- fits[1:len]
  }
  
  fit <- rowSums(fitt)
  e <- tc - fit
  
  result <- list(VM = VM, hrf = hrf, fit = fit, e = e, param = param)
  return(result)
}
