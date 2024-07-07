# Canonical DD method for Region 1

V <- nrow(dat1)

Resfit <- list()
for (i in 1:V) {
  fit <- Fit_Canonical_HRF(dat1[i,], TR, Runc, 30, 2)
  Resfit[[i]] <- fit
}

# Nonspatial Results

numstim <- 4
kbeta <- 2

# Nonspatial Beta
Nonspatial_beta <- array(0, dim = c(V, numstim * kbeta))

for (i in 1:V){
  Nonspatial_beta[i,] <- c(t(Resfit[[i]]$info$b))
}

# Fitted Amplitude, Time-to-peak, and Width
A_DD_fit <- T_DD_fit <- W_DD_fit <- array(0, dim = c(V, numstim))

for (i in 1:V){
  A_DD_fit[i, ] <- Resfit[[i]]$param[1,]
  T_DD_fit[i, ] <- Resfit[[i]]$param[2,]
  W_DD_fit[i, ] <- Resfit[[i]]$param[3,]
}

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

coords <- t(read.table("RegionData/subject-0002-run1-region1-coordinates.txt", sep = ",",
                       header = FALSE, quote = "", stringsAsFactors = FALSE))
stim_name <- c("Cue", "Expect Rating", "Stimulus", "Outcome Rating")
par(mfrow=c(2,2))

for (j in 1:4){
  for (i in (2*j-1):(2*j)){
    cols <- myColorRamp(c("white", "blue"), Nonspatial_beta[,i])
    scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                  xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
  }
}

###############################################################

# Plot Non-Spatial Amplitude

for (j in 1:4){
  cols <- myColorRamp(c("white", "red"), A_DD_fit[,j])
  scatterplot3d(coords[,1],coords[,2],coords[,3], pch = 16, color=cols,
                xlab = "x", ylab = "y", zlab = "z", box=FALSE, main = stim_name[j])
}

