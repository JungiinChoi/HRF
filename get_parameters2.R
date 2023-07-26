get_parameters2 <- function(hdrf, t) {
  # Find model parameters
  #
  # Height - h
  # Time to peak - p (in time units of TR seconds)
  # Width (at half peak) - w
  
  # Calculate Heights and Time to peak:
  n <- round(length(t) * 0.8)
  
  p <- which.max(abs(hdrf[1:n]))
  h <- hdrf[p]
  
  if (p > t[length(t)] * 0.8) {
    warning('Late time to peak')
  }

  v <- rep(0, length(hdrf))
  for (i in 1:length(hdrf)){
    if (h>0){
      v[i] <- (hdrf[i] >= (h/2))
    }else{
      v[i] <- (hdrf[i] <= (h/2))
    }
  }
  b <- which(diff(v) == -1)[1]
  v[(b + 1):length(v)] <- 0
  w <- sum(v)
  
  cnt <- p
  g <- hdrf[-1] - hdrf[-length(hdrf)]
  while (cnt > 1 && abs(g[cnt - 1]) < 0.001) {
    h <- hdrf[cnt]
    p <- cnt
    cnt <- cnt - 1
  }
  
  param <- c(h, p, w)
  return(param)
}
