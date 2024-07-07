
# Fits FIR and smooth FIR model  

# INPUTS:
# tc    - time course
# TR    - time resolution
# Runs  - experimental design
# T     - length of estimated HRF
# mode  - deterministic or stochastic
#   options:
#       0 - deterministic aproach 
#       1 - simulated annealing approach
#       Please note that when using simulated annealing approach you
#       may need to perform some tuning before use.
# 
# OUTPUTS:
# hrf   - estimated hemodynamic response function
# fit   - estimated time course
# e     - residual time course
# param - estimated amplitude, height and width

Fit_Logit2 <- function(tc, TR, Run, T, mode) {
  numstim <- length(Run)
  len <- length(Run[[1]])
  t <- seq(from = 1, to = T, by = TR)
  
  V0 <- c(1, 6, 1, 0.5, 10, 1, 15)  # initial values for logit fit: fitted values in MATLAB
  V0 <- matrix(rep(V0, numstim), nrow = numstim, byrow = TRUE)
  
  if (mode == 1 && numstim > 1) {
    cat("Multi-stimulus annealing function currently not implemented. Switching to 'deterministic mode'\n")
    mode <- 0
  }
  
  # Estimate theta (logit parameters)
  if (mode == 1) {
    cat("Stochastic Mode\n")
    
    Runs <- Run[[1]]
    # Implement Anneal_Logit function (not provided) using simulated annealing
    # [theta,HH,C,P,hrf,fit,e,param] = Anneal_Logit(V0,t,tc,Runs); 
    
  } else if (mode == 0) {
    cat("Deterministic Mode\n")
    Det_Logit_fit <- Det_Logit(V0,t,tc,Run)
    hrf = Det_Logit_fit$hrf
    fit = Det_Logit_fit$fit
    e = Det_Logit_fit$e
    param = Det_Logit_fit$param
  }
  return(list(hrf = hrf, fit = fit, e = e, param = param))
}

# Downsampling
Fit_Logit2_down <- function(tc, TR, Run, T, down) {
  numstim <- length(Run)
  len <- length(Run[[1]])
  t <- seq(from = 1, to = T, by = TR / down)
  
  V0 <- c(1, 6, 1, 0.5, 10, 1, 15)  # initial values for logit fit: fitted values in MATLAB
  V0 <- matrix(rep(V0, numstim), nrow = numstim, byrow = TRUE)
  
  Det_Logit_fit <- Det_Logit_down(V0, t, tc, Run, down)
  hrf = Det_Logit_fit$hrf
  fit = Det_Logit_fit$fit
  e = Det_Logit_fit$e
  param = Det_Logit_fit$param

  return(list(hrf = hrf, fit = fit, e = e, param = param))
}
