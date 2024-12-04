tor_make_deconv_mtx3 <- function(sf, tp, eres) {
  if (!is.list(sf)) {
    sf2 <- lapply(1:ncol(sf), function(i) sf[, i])
    sf <- sf2
  }
  
  if (length(tp) == 1) tp <- rep(tp, length(sf))
  if (length(tp) != length(sf)) stop('timepoints vectors (tp) and stick function (sf) lengths do not match!')
  
  tbefore <- if (length(args) > 2) args[[3]] else 0
  nsess <- if (length(args) > 1) args[[2]] else 1
  shiftElements <- eres
  
  if (length(args) > 4) {
    numframes <- args[[5]]
    st <- cumsum(c(1, numframes))
    en <- st[-1] - 1
    st <- st[-length(st)]
    
    DXs <- lapply(1:length(numframes), function(sess) {
      sfsess <- lapply(1:length(sf), function(i) sf[[i]][st[sess]:en[sess]])
      tor_make_deconv_mtx3(sfsess, tp, eres, args[[1]], args[[2]], 0)
    })
    
    DX <- do.call(rbind, DXs)
    DX <- DX[, -ncol(DX)]
    DX <- cbind(DX, intercept_model(rep(numframes, length(numframes))))
    
  } else {
    numtrs <- round(length(sf[[1]]) / eres)
    myzeros <- rep(0, numtrs)
    origsf <- sf
    
    sf <- lapply(sf, function(sfi) {
      Snumtrs <- length(sfi) / eres
      if (Snumtrs != round(Snumtrs)) warning(paste('sf{', which(sf == sfi), '}: length not evenly divisible by eres.'))
      if (numtrs != Snumtrs) warning(paste('sf{', which(sf == sfi), '}: different length than sf{1}.'))
      inums <- ceiling(which(sfi > 0) / eres)
      inums[inums == 0] <- 1
      new_sf <- myzeros
      new_sf[inums] <- 1
      new_sf
    })
    
    DX <- matrix(0, nrow = numtrs, ncol = sum(tp) + if (nsess < 2) 1 else nsess)
    index <- 1
    
    for (i in 1:length(sf)) {
      if (tbefore != 0) {
        for (j in tbefore:1) {
          mysf <- c(sf[[i]][(j + 1):numtrs], rep(0, j))
          DX[, index] <- mysf
          index <- index + 1
        }
      }
      
      DX[, index] <- sf[[i]]
      index <- index + 1
      inums <- which(sf[[i]] == 1)
      
      for (j in 2:tp[i]) {
        inums <- inums + 1
        reg <- myzeros
        reg[inums] <- 1
        reg <- head(reg, numtrs)
        while (length(reg) < nrow(DX)) reg <- c(reg, 0)
        DX[, index] <- reg
        index <- index + 1
      }
    }
    
    if (nsess < 2) {
      DX[, ncol(DX)] <- 1
    } else {
      scanlen <- nrow(DX) / nsess
      if (round(scanlen) != scanlen) warning('Model length is not an even multiple of scan length.')
      X <- matrix(0, nrow(DX), nsess)
      for (startimg in seq(1, nrow(DX), by = scanlen)) {
        X[startimg:(startimg + scanlen - 1), (startimg - 1) %/% scanlen + 1] <- 1
      }
      DX <- cbind(DX, X)
    }
  }
  
  list(DX = DX, sf = sf)
}
