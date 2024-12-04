# Load required libraries
library(ggplot2)

# Define the regions and initialize empty lists
regions <- 140:147
hrf_lists <- lapply(regions, function(x) list(hrf = numeric(), mean = NULL, sd = NULL))

# Combine H matrix
H <- cbind(CanonicalBasisSet(TR)$h, CanonicalBasisSet(TR)$dh)

# Calculate HRF for each region
for (n in regions) {
  region_index <- which(regions == n)
  for (i in 1:reglen[n]) {
    b2 <- matrix(Nonspatial_beta_spline[n, i, ], nrow = numstim, ncol = 5, byrow = TRUE)
    b2[,1:2] <- 4*b2[,1:2]
    hrf_tmp <- B %*% t(matrix(b2[3,], nrow = 1))
    hrf_lists[[region_index]]$hrf <- c(hrf_lists[[region_index]]$hrf, c(hrf_tmp))
  }
}

# Compute means and standard deviations
for (i in seq_along(hrf_lists)) {
  hrf_matrix <- matrix(hrf_lists[[i]]$hrf, nrow = reglen[regions[i]], ncol = 64, byrow = TRUE)
  hrf_lists[[i]]$mean <- apply(hrf_matrix, 2, mean)
  hrf_lists[[i]]$sd <- apply(hrf_matrix, 2, sd)
}

# Prepare data frame for plotting
DD_spat_df <- data.frame(
  region = factor(rep(seq_along(regions), each = 64)),
  mean_value = unlist(lapply(hrf_lists, `[[`, "mean")),
  sd = unlist(lapply(hrf_lists, `[[`, "sd")),
  t = rep(1:64, length(regions))
)



# Create the plot
ggplot(DD_spat_df, aes(x = t, y = mean_value, color = region)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_value - 0.5 * sd, ymax = mean_value + 0.5 * sd, fill = region), alpha = 0.2, color = NA) +  
  labs(title = "Spatial HRF: Region 175",
       x = "Time",
       y = "y",
       fill = "Subject",
       color = "Subject") +
  theme_minimal() +
  theme(legend.position = "right")


for (i in 1:512){
  if (DD_spat_df$t[i] >= 10 & DD_spat_df$t[i] >= 30){
    DD_spat_df$mean_value[i] <- 2* DD_spat_df$mean_value[i]
  }
}
