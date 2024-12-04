# Visualize

library(rgl)

myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

myColorRamp_minmax <- function(colors, values, min, max) {
  v <- (values - min)/(max-min)
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

Nonspatial_beta_DD_merge <- matrix(ncol = 8)
coords_merge <- matrix(ncol = 3)

for (i in 1:numreg){
  Nonspatial_beta_DD_merge <- rbind(Nonspatial_beta_DD_merge, 
                                    Nonspatial_beta_DD[i,1:reglen[i],])
  coords_merge <- rbind(coords_merge, coord_list[[i]])
}

Nonspatial_beta_DD_merge<- Nonspatial_beta_DD_merge[-1,]
coords_merge <- coords_merge[-1,]

stim_name <- c("Cue", "Expect Rating", "Stimulus", "Outcome Rating")
par(mfrow=c(2,2))
for (j in 1:4){
  for (i in (2*j-1):(2*j)){
    cols <- myColorRamp(c("white", "blue"), Nonspatial_beta_DD_merge[,i])
    scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
  }
}

###############################################################

# Plot Non-Spatial Amplitude

A_DD_fit_merge <- matrix(ncol = 4)

for (i in 1:numreg){
  A_DD_fit_merge <- rbind(A_DD_fit_merge, 
                          A_DD_fit[i,1:reglen[i],])
}

A_DD_fit_merge<- A_DD_fit_merge[-1,]

for (j in 1:4){
  cols <- myColorRamp(c("white", "red"), A_DD_fit_merge[,j])
  scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

#############################

# Visualize

spatial_beta_DD_merge <- matrix(ncol = 8)

for (i in 1:numreg){
  spatial_beta_DD_merge <- rbind(spatial_beta_DD_merge, 
                                 beta_spatial[[i]])
}

spatial_beta_DD_merge <- spatial_beta_DD_merge[-1,]


stim_name <- c("Cue", "Expect Rating", "Stimulus", "Outcome Rating")
par(mfrow=c(2,2))
for (j in 1:4){
  for (i in (2*j-1):(2*j)){
    cols <- myColorRamp(c("white", "blue"), spatial_beta_DD_merge[,i])
    scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
  }
}

###############################################################


# Plot Spatial Amplitude

A_DD_spatial <- array(0, dim = c(numreg, max(reglen), numstim))

TR <- 0.460
H <- cbind(CanonicalBasisSet(TR)$h, CanonicalBasisSet(TR)$dh)

for (i in 1:numreg){
  for (k in 1:reglen[i]){
    b <- matrix(beta_spatial[[i]][k,], nrow = numstim, ncol = 2, byrow = TRUE)
    for (j in 1:numstim) {
      hrf_tmp <- H %*% t(matrix(b[j, ], nrow = 1))
      A_DD_spatial[i,k,j] <- get_parameters2(hrf_tmp, 1:64)[1]
    }
  }
}

# Plot Spatial Amplitude

A_DD_spatial_merge <- matrix(ncol = 4)

for (i in 1:numreg){
  A_DD_spatial_merge <- rbind(A_DD_spatial_merge, 
                              A_DD_spatial[i,1:reglen[i],])
}

A_DD_spatial_merge<- A_DD_spatial_merge[-1,]

for (j in 1:4){
  cols <- myColorRamp(c("white", "red"), A_DD_spatial_merge[,j])
  scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

###############
#Spline

Nonspatial_beta_spline_merge <- matrix(ncol = 20)

for (i in 1:numreg){
  Nonspatial_beta_spline_merge <- rbind(Nonspatial_beta_spline_merge, 
                                    Nonspatial_beta_spline[i,1:reglen[i],])
}

Nonspatial_beta_spline_merge <- Nonspatial_beta_spline_merge[-1,]

stim_name <- c("Cue", "Expect Rating", "Stimulus", "Outcome Rating")
par(mfrow=c(1,5))
for (j in 1:4){
  for (i in (5*j-4):(5*j)){
    cols <- myColorRamp(c("white", "blue"), Nonspatial_beta_spline_merge[,i])
    scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
  }
}


#Spline

spatial_beta_spline_merge <- matrix(ncol = 20)

for (i in 1:numreg){
  spatial_beta_spline_merge <- rbind(spatial_beta_spline_merge, 
                                        beta_spatial_spline[[i]])
}

spatial_beta_spline_merge <- spatial_beta_spline_merge[-1,]

stim_name <- c("Cue", "Expect Rating", "Stimulus", "Outcome Rating")
par(mfrow=c(1,5))
for (j in 1:4){
  for (i in (5*j-4):(5*j)){
    cols <- myColorRamp(c("white", "blue"), spatial_beta_spline_merge[,i])
    scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
  }
}

###############################################################

# Plot Non-Spatial Amplitude

A_spline_fit_merge <- matrix(ncol = 4)

for (i in 1:numreg){
  A_spline_fit_merge <- rbind(A_spline_fit_merge, 
                          A_spline_fit[i,1:reglen[i],])
}

A_spline_fit_merge<- A_spline_fit_merge[-1,]

for (j in 1:4){
  cols <- myColorRamp(c("white", "red"), A_spline_fit_merge[,j])
  scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

# Plot Spatial Amplitude

A_spline_spatial <- array(0, dim = c(numreg, max(reglen), numstim))

TR <- 0.460
H <- cbind(CanonicalBasisSet(TR)$h, CanonicalBasisSet(TR)$dh)

for (i in 1:numreg){
  for (k in 1:reglen[i]){
    b <- matrix(beta_spatial_spline[[i]][k,], nrow = numstim, ncol = 5, byrow = TRUE)
    for (j in 1:numstim) {
      hrf_tmp <- H %*% t(matrix(b[j, ], nrow = 1))
      A_spline_spatial[i,k,j] <- get_parameters2(hrf_tmp, 1:64)[1]
    }
  }
}

# Plot Spatial Amplitude

A_spline_spatial_merge <- matrix(ncol = 4)

for (i in 1:numreg){
  A_DD_spatial_merge <- rbind(A_DD_spatial_merge, 
                              A_DD_spatial[i,1:reglen[i],])
}

A_DD_spatial_merge<- A_DD_spatial_merge[-1,]

for (j in 1:4){
  cols <- myColorRamp(c("white", "red"), A_DD_spatial_merge[,j])
  scatterplot3d(coords_merge[,1],coords_merge[,2],coords_merge[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

