load("/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR/DD_1.RData")

stim_name <- c("Cue", "LF", "RF", "LH", "RH", "Tongue")


A_DD_fit <- A_DD_fit_1
A_DD_mean <- apply(A_DD_fit, c(2,3), mean)
A_DD_spat_mean <- cbind(A_DD_spat_mean, coords)
write.csv(A_DD_mean, file = "/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR/Spatial_DD_A.csv", 
          row.names = FALSE)


nonspatial_DD_mean <- apply(Nonspatial_beta_1, c(2,3), mean)

write.csv(A_DD_spat_mean, file = "/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR/Nonspatial_DD_A.csv", 
          row.names = FALSE)

write.csv(A_DD_spat_mean, file = "/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR/spatial_DD_A.csv", 
          row.names = FALSE)

write.csv(A_spline_mean, file = "/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR/Nonspatial_sp_A.csv", 
          row.names = FALSE)

A_spline_mean <- apply(A_spline_spatial, c(2,3), mean)
write.csv(A_spline_mean, file = "/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR/spatial_sp3_A.csv", 
          row.names = FALSE)

write.csv(A_sFIR_fit, file = "/Users/user/Documents/JHU/research/HRF_Est/HCM_MOTOR/Nonspatial_sf_A.csv", 
          row.names = FALSE)


A_sp_spat <- matrix(0, nrow = sum(reglen), ncol = 7)
ind <- 1
for (i in 1:numreg){
  for (j in 1:reglen[i]){
    A_sp_spat[ind,1:4] <- A_sp_spatial[i,j,]
    A_sp_spat[ind,5:7] <- coord_list[[i]][j,]
    ind <- ind + 1
  }
}
write.csv(e_DD_nonspat, file = "/Users/user/Documents/JHU/research/HRF_Est/Pain/e_DD_nonspat.csv", 
          row.names = FALSE)
