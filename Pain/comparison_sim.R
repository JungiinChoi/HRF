# Simulation to compare Gaussian Kernel Smoothing vs TPRS vs Shrinkage methods

# Load Rscripts

source("../Fit_Canonical_HRF.R")
source("../Fit_Spline.R")
source("../get_parameters2.R")
source("../spm_hrf.R")
source("../Det_Logit.R")
source("../tor_make_deconv_mtx3.R")

# Load necessary libraries
library(ggplot2)
library(tidyverse)
library(sp)
library(gstat)
library(raster)
library(stats)
library(mgcv)  # for thin plate regression spline
library(akima)
library(Matrix)

# Model: Canonical DD (two betas)

# I. 2D Square Region

# Define the test function f(x, z)
f_function_1 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (3.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.2)^2 / sigma_x^2) - ((z - 0.3)^2 / sigma_z^2))
  term2 <- (2.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.7)^2 / sigma_x^2) - ((z - 0.8)^2 / sigma_z^2))
  return(term1 + term2)
}

f_function_2 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (0.3 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.1)^2 / sigma_z^2) - ((z - 0.8)^2 / sigma_x^2))
  term2 <- (0.6 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.9)^2 / sigma_z^2) - ((z - 0.2)^2 / sigma_x^2))
  return(term1 + term2)
}

# Set parameters
sigma_x <- 0.3
sigma_z <- 0.4
n_points <- 500  # number of points to generate

# Generate random (x, z) points within the unit square
set.seed(123)
grid <- expand.grid(x = seq(0, 1, length.out = sqrt(n_points)), y = seq(0, 1, length.out = sqrt(n_points)))

x <- grid[,1]
z <- grid[,2]

# Create a data frame from the vectors
data <- data.frame(x1 = x, x2 = z)

# Use ggplot to plot the points
ggplot(data, aes(x = x1, y = x2)) +
  geom_point(color = "blue") +  # Plot points
  labs(title = "Voxels", x = "X1", y = "X2") +  # Add labels
  theme_minimal()  # Use a clean minimal theme

beta1 <- f_function_1(x, z, sigma_x, sigma_z)
beta2 <- f_function_2(x, z, sigma_x, sigma_z)

# Create a grid for plotting the true function
grid_size <- 50
x_grid <- seq(0, 1, length.out = grid_size)
z_grid <- seq(0, 1, length.out = grid_size)
beta1_grid <- outer(x_grid, z_grid, function(x, z) f_function_1(x, z, sigma_x, sigma_z))
beta2_grid <- outer(x_grid, z_grid, function(x, z) f_function_2(x, z, sigma_x, sigma_z))

# Draw contour plots
par(mfrow = c(1, 2))  # Set up the plotting window for two plots side by side

# Contour plot of the true function
contour(x_grid, z_grid, beta1_grid, main = "True Beta1", xlab = "x1", ylab = "x2")
contour(x_grid, z_grid, beta2_grid, main = "True Beta2", xlab = "x1", ylab = "x2")

# Case I: Large noise (fmri data has large noise!)
# 1. Independent(Non-spatial) OLS fitting

# Settings
TR = 1
T = 40
t = seq(1, T, by = TR)    # samples at which to get Logit HRF Estimate
FWHM = 4                  # FWHM for residual scan
pval = 0.01
df = 600
alpha = 0.001
len = 640

# function generating error 
noise_arp <- function(n, phi = c(0.5, 0.1), sigma = 1, e = rnorm(n, 0, sigma)){
  p <- length(phi)
  w <- numeric(n)
  
  for (i in 1:n) {
    if (i <= p) {
      w[i] <- e[i]
    } else {
      arpart <- sum(phi * w[(i-1):(i-p)], na.rm = TRUE)
      w[i] <- arpart + e[i]
    }
  }
  return(w)
}

# Plot the hrf
beta_true_1 <- seq(0.317,9.044,length = 30)
beta_true_2 <- seq(0.0104,1.53,length = 30)
hrfdf <- data.frame()
hrfmat <- matrix(nrow = length(beta_true_1), ncol = 40)
for (i in 1:length(beta_true_1)){
  true_hrf <- CanonicalBasisSet(TR)$h * beta_true_1[i] + CanonicalBasisSet(TR)$dh * beta_true_2[i]
  hrfdf <- rbind(hrfdf, data.frame(index = rep(i,40), x = 0:39, y = true_hrf))
}


hrfdf$index <- factor(hrfdf$index)

ggplot(hrfdf, aes(x, y, group = index, colour = index)) +
  labs(title = "HRF for Different Betas", x = "t", y = "y") +
  geom_line() 

# Onset
b <- 1
R <- c(13, 14, 29, 44, 107, 125, 160, 171, 174, 190, 191, 
       206, 215, 232, 237, 262, 277, 292, 296, 346, 354, 
       367, 375, 382, 384, 398, 409, 462, 469, 475, 501, 
       520, 527, 566, 577, 629)
Run <- rep(0, 640)
Run[R] <- 1
Runc <- list(Run)
onset <- rep(-100,640)
onset[R] <- R

# Simulations
set.seed(906)
tc_mg <- c()

sim_k <- 5
error_sig <- 5
beta_DD <- array(0, dim = c(length(beta1), sim_k, 2))

for (k in 1:length(beta1)){
  true_sig <-  beta1[k] * conv(Run, CanonicalBasisSet(TR)$h)[1:len] + beta2[k] * conv(Run, CanonicalBasisSet(TR)$dh)[1:len] 
  
  xsecs <- 0:40
  tc_mat <- matrix(0,nrow = len, ncol = 5)
  
  for (i in 1:sim_k){
    tc_noise <- noise_arp(n = len, phi = c(0.3, 0))
    tc <- true_sig + error_sig * tc_noise
    if (i ==1){
      tc_mg <- c(tc_mg, tc)
    }
    tc_mat[,i] <- tc
    
    #Canonical HRF
    Canonical_HRF_fitted <- Fit_Canonical_HRF(tc, TR, Runc, T, 2)
    beta_DD[k,i,] <- Canonical_HRF_fitted$info$b[1,]
  }
  print(k)
}

# Visualize timecourse
HRF_df <- data.frame(t = 1:len, tc = tc_mg[1:640], onset = onset)
#HRF_df$index = factor(HRF_df$index)

ggplot(HRF_df) +
  geom_line(aes(x=t, y=tc)) +
  theme_minimal() +
  labs(title = "Convoluted time course", x = "", y = "") +
  geom_segment(aes(x = onset, y = -20, xend = onset, yend = -15)) +
  coord_cartesian(xlim = c(0, 640))



# Estimate beta independently

ind_beta1 <- apply(beta_DD[,,1],1,mean)
ind_beta2 <- apply(beta_DD[,,2],1,mean)
ind_dist <- (beta1 - ind_beta1)^2 + 
  (beta2 - ind_beta2)^2

df_fitted <- data.frame(x1 = x, x2 = z, ind_beta1 = ind_beta1, ind_beta2 = ind_beta2)

#create plot

fld <- with(df_ind, interp(x = x1, y = x2, z = ind_beta1))
fld$x <- c(fld$x, 2.000001)
fld$y <- c(fld$y, 2.000001)
fld$z <- rbind(cbind(fld$z, rep(10, 40)), rep(10,41))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "Independent Beta 1",
               key.title = title(main = "beta1 ", cex.main = 1))

fld <- with(df_ind, interp(x = x1, y = x2, z = ind_beta2))
fld$x <- c(fld$x, 1.000001)
fld$y <- c(fld$y, 1.000001)
fld$z <- rbind(cbind(fld$z, rep(2, 40)), rep(2,41))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "Independent Beta 2",
               key.title = title(main = "beta2", cex.main = 1))

contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(ind_beta1, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
        main = "True Beta1", xlab = "x1", ylab = "x2")
contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(ind_beta2, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
        main = "True Beta2", xlab = "x1", ylab = "x2")

fld <- with(df_ind, interp(x = x1, y = x2, z = ind_dist))
fld$x <- c(fld$x, 1.000001)
fld$y <- c(fld$y, 1.000001)
fld$z <- rbind(cbind(fld$z, rep(1.2, 40)), rep(1.2,41))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "L2 Distance From True Beta",
               key.title = title(main = "d", cex.main = 1))


# 2. TPRS (Nonparametric)

beta_TPRS <- list()

k <- 1
V <- length(x)
for (k_gam in c(20,90)){
  gamma_fitted <- matrix(0, nrow = V, ncol = 2 * k_gam)
  
  gam_df <- data.frame(beta = ind_beta1, x = x, 
                       y = z)
  gamfit <- gam(beta ~ s(x,y, k = k_gam), data = gam_df, method = "REML")
  Blist <- as.matrix(model.matrix(gamfit))
  
  X <- Canonical_HRF_fitted$info$X[,-1]
  X_ext <- bdiag(replicate(V, X, simplify = FALSE))
  
  i_perm <- 1:(2 * V)
  j_perm <- c()
  for (i in 1:V){
    j_perm <- c(j_perm, (0:(2-1))*V+i)
  }
  perm <- sparseMatrix(i = i_perm, j = j_perm, x = 1) 
  
  B <- bdiag(replicate(2, Blist, simplify = FALSE))
  PB <- perm %*% B
  y_ext <- tc_mg
  
  XtX <- t(X_ext) %*% X_ext
  BPXtX <- t(PB) %*% XtX
  A <- BPXtX %*% PB
  Xty_ext <- t(X_ext) %*% y_ext
  b <- t(PB) %*% Xty_ext
  
  gamma_fitted[k,] <- Matrix::solve(A, b)[,1]
  beta_TPRS[[k]] <- matrix(0, V, 2)
  for (i in 1:(2)){
    beta_TPRS[[k]][,i] <- (Blist %*% gamma_fitted[k,(1+(i-1)*k_gam):(i*k_gam)])[,1]
  }
  print(k)
  k <- k+1
}


ind <- 1
tprs_beta1 <- beta_TPRS[[ind]][,1]
tprs_beta2 <- beta_TPRS[[ind]][,2]
tprs_dist <- (beta1 - tprs_beta1)^2 + 
  (beta2 - tprs_beta2)^2

#create plot
df_fitted$tprs_beta1 <- beta_TPRS[[ind]][,1]
df_fitted$tprs_beta2 <- beta_TPRS[[ind]][,2]

fld <- with(df_fitted, interp(x = x1, y = x2, z = tprs_beta1))
fld$x <- c(fld$x, 1.000001)
fld$y <- c(fld$y, 1.000001)
fld$z <- rbind(cbind(fld$z, rep(10, 40)), rep(10,41))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "TPRS Beta 1",
               key.title = title(main = "beta1 ", cex.main = 1))

fld <- with(df_fitted, interp(x = x1, y = x2, z = tprs_beta2))
fld$x <- c(fld$x, 1.000001)
fld$y <- c(fld$y, 1.000001)
fld$z <- rbind(cbind(fld$z, rep(2, 40)), rep(2,41))

filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "TPRS Beta 2",
               key.title = title(main = "beta2", cex.main = 1))

contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(tprs_beta1, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
        main = "Kernel Smoothing Beta1", xlab = "x1", ylab = "x2")
contour(seq(0, 1, length.out = sqrt(n_points)), seq(0, 1, length.out = sqrt(n_points)), 
        matrix(tprs_beta2, nrow = length(seq(0, 1, length.out = sqrt(n_points))), byrow = T), 
        main = "Kernel Smoothing Beta2", xlab = "x1", ylab = "x2")

fld <- with(df_ind, interp(x = x1, y = x2, z = tprs_dist))
fld$x <- c(fld$x, 1.000001)
fld$y <- c(fld$y, 1.000001)
fld$z <- rbind(cbind(fld$z, rep(1.2, 40)), rep(1.2,41))
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "X1",
               ylab = "X2",
               main = "L2 Distance",
               key.title = title(main = "d", cex.main = 1))

ggplot(data.frame(k = rep(k_vec,2), ll = c(ll,bic), criteria = as.factor(c(rep("aic", 15), rep("bic", 15)))), 
       aes(x=k, y=ll, group = criteria)) +
  geom_line(aes(color = criteria))+
  geom_point(aes(color = criteria)) +
  theme_minimal() +
  labs(title = "IC for different number of spline basis",
       x = "number of basis",
       y = "Information")

# 3. Gaussian Kernel Smoothing 



# 4. Shrinkage method




# II. 2D Multiple Parcells


# Model: Canonical DD (two betas)

# I. 2D Square Region

# Define the test function f(x, z)
f_func_11 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (3.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.2)^2 / sigma_x^2) - ((z - 0.3)^2 / sigma_z^2))
  term2 <- (2.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.7)^2 / sigma_x^2) - ((z - 0.8)^2 / sigma_z^2))
  return(term1 + term2)
}

f_func_12 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (0.3 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.1)^2 / sigma_z^2) - ((z - 0.8)^2 / sigma_x^2))
  term2 <- (0.6 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.9)^2 / sigma_z^2) - ((z - 0.2)^2 / sigma_x^2))
  return(term1 + term2)
}

f_func_21 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (4 / (pi * sigma_x * sigma_z * 3)) * exp(-((x - 0.2)^2 / 2 / sigma_x^2) - ((z - 0.3)^2 / 3 / sigma_z^2))
  term2 <- (1 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.7)^2 / sigma_x^2) - ((z - 0.8)^2 / sigma_z^2))
  return(term1 + term2)
}

f_func_22 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (0.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.1)^2 / sigma_z^2) - ((z - 0.8)^2 / sigma_x^2))
  return(term1)
}

f_func_31 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (1 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.2)^2 / 2/ sigma_x^2) - ((z - 0.3)^2 / 6 / sigma_z^2))
  term2 <- (2.5 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.7)^2 / 5 / sigma_x^2) - ((z - 0.8)^2 / 2/ sigma_z^2))
  return(term1 + term2)
}

f_func_32 <- function(x, z, sigma_x, sigma_z) {
  term2 <- (0.6 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.9)^2 / 3 / sigma_z^2) - ((z - 0.2)^2 / 2 / sigma_x^2))
  return(term2)
}

f_func_41 <- function(x, z, sigma_x, sigma_z) {
  term1 <- (4 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.2)^2 / 6 / sigma_x^2) - ((z - 0.3)^2 / 1 / sigma_z^2))
  return(term1)
}

f_func_42 <- function(x, z, sigma_x, sigma_z) {
  term2 <- (0.8 / (pi * sigma_x * sigma_z)) * exp(-((x - 0.9)^2 / sigma_z^2) - ((z - 0.2)^2 / sigma_x^2))
  return(term2)
}

f_func_1 <- function(x, z, sigma_x, sigma_z){
  res <- rep(0,length(x))
  for (i in 1:length(x)){
    if (x[i]<1 & z[i]<1){
      res[i] <- f_func_11(x[i], z[i], sigma_x, sigma_z)
    } else if (x[i]<1){
      res[i] <- f_func_21(x[i], z[i]-1, sigma_x, sigma_z)
    } else if (z[i]<1){
      res[i] <- f_func_31(x[i]-1, z[i], sigma_x, sigma_z)
    }else{
      res[i] <- f_func_41(x[i]-1, z[i]-1, sigma_x, sigma_z)
    }
  }
  return(res)
}
f_func_2 <- function(x, z, sigma_x, sigma_z){
  res <- rep(0,length(x))
  for (i in 1:length(x)){
    if (x[i]<1 & z[i]<1){
      res[i] <- f_func_12(x[i], z[i], sigma_x, sigma_z)
    } else if (x[i]<1){
      res[i] <- f_func_22(x[i], z[i]-1, sigma_x, sigma_z)
    } else if (z[i]<1){
      res[i] <- f_func_32(x[i]-1, z[i], sigma_x, sigma_z)
    }else{
      res[i] <- f_func_42(x[i]-1, z[i]-1, sigma_x, sigma_z)
    }
  }
  return(res)
}

# Set parameters
sigma_x <- 0.3
sigma_z <- 0.4
n_points <- 1000  # number of points to generate

# Generate random (x, z) points within the unit square
set.seed(123)
grid <- expand.grid(x = seq(0, 2, length.out = sqrt(n_points)), y = seq(0, 2, length.out = sqrt(n_points)))

x <- grid[,1]
z <- grid[,2]

beta1 <- f_func_1(x, z, sigma_x, sigma_z)
beta2 <- f_func_2(x, z, sigma_x, sigma_z)

# Create a grid for plotting the true function
grid_size <- 100
x_grid <- seq(0, 2, length.out = grid_size)
z_grid <- seq(0, 2, length.out = grid_size)
beta1_grid <- outer(x_grid, z_grid, function(x, z) f_func_1(x, z, sigma_x, sigma_z))
beta2_grid <- outer(x_grid, z_grid, function(x, z) f_func_2(x, z, sigma_x, sigma_z))

# Draw contour plots
par(mfrow = c(1, 2))  # Set up the plotting window for two plots side by side

# Contour plot of the true function
contour(x_grid, z_grid, beta1_grid, main = "True Beta1", xlab = "x1", ylab = "x2")
contour(x_grid, z_grid, beta2_grid, main = "True Beta2", xlab = "x1", ylab = "x2")


# III. Irregular surface
# Load necessary libraries
library(sp)
library(ggplot2)

# Define the boundary points of the shape based on the visual of the image you provided
boundary_points <- data.frame(
  x = c(0.1, 0.35, 0.4, 0.6, 0.6, 0.75, 0.85, 0.75, 0.6, 0.4, 0.3, 0.15),
  y = c(0.5, 0.45, 0.85, 0.75, 0.6, 0.4, 0.2, 0.1, 0.2, 0.15, 0.3, 0.4)
)

# Create a Polygon object
polygon_region <- Polygon(boundary_points)

# Create a function to check whether a point (x, y) is inside the polygon
is_inside_region <- function(x, y) {
  # Create SpatialPolygons object
  sp_polygon <- SpatialPolygons(list(Polygons(list(polygon_region), ID = 1)))
  
  # Check if the point (x, y) is inside the polygon
  point <- SpatialPoints(data.frame(x = x, y = y))
  inside <- over(point, sp_polygon)
  
  return(!is.na(inside))  # Return 1 if inside, 0 if outside
}

# Test the function with some points
x_test <- c(0.2, 0.5, 0.9)
y_test <- c(0.3, 0.8, 0.5)
result <- sapply(1:length(x_test), function(i) is_inside_region(x_test[i], y_test[i]))
print(result)  # 1 if inside, 0 if outside

# Visualization of the polygon and grid of points
x_vals <- seq(0, 1, length.out = 50)
y_vals <- seq(0, 1, length.out = 50)
grid_data <- expand.grid(x = x_vals, y = y_vals)

# Apply the function to check if points are inside the region
grid_data$z <- apply(grid_data, 1, function(row) is_inside_region(row['x'], row['y']))

# Plot the shape and the points
ggplot() +
  geom_polygon(data = boundary_points, aes(x = x, y = y), fill = "lightblue", color = "black") +
  geom_point(data = grid_data, aes(x = x, y = y, color = factor(z)), size = 0.5) +
  scale_color_manual(values = c("white", "blue"), labels = c("Outside", "Inside")) +
  labs(title = "Irregular Shape with Inside Region", color = "Region") +
  theme_minimal() +
  coord_fixed(ratio = 1)

# Gaussian Kernel smoothing 

# Load necessary libraries
library(MASS)  # for gaussian kernel smoothing
library(ggplot2)

# Function to calculate the Gaussian kernel weights
gaussian_kernel <- function(x1, y1, x2, y2, sigma) {
  dist_sq <- (x1 - x2)^2 + (y1 - y2)^2
  return(exp(-dist_sq / (2 * sigma^2)))
}

# Function to perform Gaussian kernel smoothing
gaussian_kernel_smoothing <- function(coords, values, sigma) {
  smoothed_values <- numeric(nrow(coords))
  
  for (i in 1:nrow(coords)) {
    weights <- numeric(nrow(coords))
    
    # Calculate weights for all points based on their distance to the i-th point
    for (j in 1:nrow(coords)) {
      weights[j] <- gaussian_kernel(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2], sigma)
    }
    
    # Weighted average of the values
    smoothed_values[i] <- sum(weights * values) / sum(weights)
  }
  
  return(smoothed_values)
}

# Example: Simulated 2D data with coordinates and values
set.seed(42)
n <- 200  # number of points
coords <- data.frame(x =grid_data$x, y = grid_data$y)  # random 2D coordinates
values <- grid_data$z

# Apply Gaussian kernel smoothing with a given sigma
sigma <- 0.01  # smoothing parameter
smoothed_values <- gaussian_kernel_smoothing(coords, values, sigma)

# Visualize the original and smoothed values
ggplot() +
  geom_point(aes(x = coords$x, y = coords$y, color = values), size = 2) +
  ggtitle("Original Data") +
  scale_color_manual(values = c("white", "blue"))+
  theme_minimal()

ggplot() +
  geom_point(aes(x = coords$x, y = coords$y, color = smoothed_values), size = 2) +
  ggtitle("Smoothed Data") +
  scale_color_gradient(low = "white", high = "blue") +
  theme_minimal()

plot(1:529, Blist[,2], type="l")
