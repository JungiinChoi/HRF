# INPUTS:

# V0 - Initial values for IL-model
# tc - time course data
# TR - time resolution
# Runc - a cell array containing stick functions (one for each condition)
# down - downsampling factor

# OUTPUTS:

# VM - final parameter values for IL-model
# hrf - estimated HRFs for each condition
# fit - estimated time course

Fit_Logit_allstim_2 <- function(V0, tc, TR, Runc, down) {
  # Initial values
  numstim <- length(Runc)
  len <- length(Runc[[1]])
  
  LB <- c(0.05, 1, 0, 0.05, 5, 0.05, 8)  # Previous Lower bounds for parameters
  UB <- c(2, 10, 2, 2, 15, 1, 20)
  
  LB <- rep(LB, numstim)
  UB <- rep(UB, numstim)
  V0 <- rep(V0, numstim)
  
  # Find optimal values
  result <- optim(par = V0, fn = cost_allstim_2, lower = LB, upper = UB, 
                  TR = TR, tc = tc, Runc = Runc, down = down, 
                  control = list(maxit = 100000, trace = 0))
  VM <- result$par
  
  # Use optimal values to fit hemodynamic response functions
  t <- seq(0, 30, by = 1 / down)
  hrf <- matrix(0, nrow = length(t), ncol = numstim)
  fitt <- matrix(0, nrow = len, ncol = numstim)
  
  for (g in 1:numstim) {
    hrf[, g] <- Get_Logit(VM[(g - 1) * 7 + 1:(g * 7)], t)  # Calculate HRF estimate (fit, given theta)
  }
  
  for (g in 1:numstim) {
    fits <- conv(Runc[[g]], hrf[, g], type = "open")
    fitt[, g] <- fits[1:len]
  }
  
  fit <- rowSums(fitt)
  
  return(list(VM = VM, hrf = hrf, fit = fit))
}

# You would need to define the following functions:
cost_allstim_2 <- function(par, TR, tc, Runc, down) {
  # Define your cost function here
}

Get_Logit <- function(theta, t) {
  # Define the Get_Logit function here
}

