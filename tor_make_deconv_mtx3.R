tor_make_deconv_mtx3 <- function(sf, tp, eres) {
  docenter <- 0
  
  if (!is.list(sf)) {
    sf <- lapply(1:ncol(sf), function(i) sf[, i])
  }
  
  if (length(tp) == 1) {
    tp <- rep(tp, length(sf))
  }
  if (length(tp) != length(sf)) {
    stop("timepoints vectors (tp) and stick function (sf) lengths do not match!")
  }
  
  tbefore <- 0
  nsess <- length(sf)
  
  shiftElements <- eres
  
  # downsample sf to number of TRs
  numtrs <- round(length(sf[[1]]) / eres)
  myzeros <- numeric(numtrs)
  origsf <- sf
  
  for (i in seq_along(sf)) {
    Snumtrs <- length(sf[[i]]) / eres
    if (Snumtrs != round(Snumtrs)) {
      warning(paste("sf[", i, "]: length not evenly divisible by eres."))
    }
    if (numtrs != Snumtrs) {
      warning(paste("sf[", i, "]: different length than sf[[1]]."))
    }
    
    inums <- which(sf[[i]] > 0) / eres  # convert to TRs
    inums <- ceiling(inums)  # nearest TR; use ceiling to avoid issues with onsets between TRs
    inums[inums == 0] <- 1  # never use 0th element
    sf[[i]] <- myzeros
    sf[[i]][inums] <- 1  # always use 1 for sf
  }
  
  # make deconvolution matrix DX
  index <- 1
  DX <- matrix(0, nrow = numtrs, ncol = length(tp) + sum(tp))
  for (i in seq_along(sf)) {
    if (tbefore != 0) {
      for (j in tbefore:1) {
        mysf <- c(sf[[i]][j + 1:length(sf[[i]])], rep(0, j))
        DX[, index] <- mysf
        index <- index + 1
      }
    }
    
    DX[, index] <- sf[[i]]
    index <- index + 1
    inums <- which(sf[[i]] == 1)
    
    for (j in 2:tp[i]) {
      inums <- inums + 1  # + 1 because we've downsampled already.  + shiftElements;
      reg <- myzeros
      reg[inums] <- 1
      reg <- reg[1:numtrs]
      while (length(reg) < nrow(DX)) {
        reg <- c(reg, 0)  # add 0's if too short
      }
      DX[, index] <- reg
      index <- index + 1
    }
  }
  
  # add intercept
  if (nsess < 2) {
    DX <- cbind(DX, rep(1, numtrs))
  } else {
    index <- 1
    scanlen <- nrow(DX) / nsess
    if (round(scanlen) != scanlen) {
      warning("Model length is not an even multiple of scan length.")
    }
    X <- matrix(0, nrow = nrow(DX), ncol = nsess)
    for (startimg in seq(1, nrow(DX), scanlen)) {
      X[startimg:(startimg + scanlen - 1), index] <- 1
      index <- index + 1
    }
    DX <- cbind(DX, X)
  }
  
  if (docenter) {
    # center columns (not intercepts)
    wh <- 1:(ncol(DX) - nsess)
    DX[, wh] <- DX[, wh] - matrix(colMeans(DX[, wh]), nrow(DX), length(wh), byrow = TRUE)
  }
  
  return(list(DX = DX, sf = sf))
}
