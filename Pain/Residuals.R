# Gaussian Kernel Smoothing

# Load necessary library
library(stats)  # For ksmooth

# Example data setup
set.seed(123)
data_length <- 300
x <- seq(1, data_length)  # Example x-coordinates
y <- seq(1, data_length)  # Example y-coordinates
z <- rnorm(data_length)   # Example 3D vector values (dependent variable)

# Perform Gaussian kernel smoothing
# For demonstration, we'll use ksmooth with a Gaussian kernel on the x-coordinates
smoothed_data <- ksmooth(x, z, kernel = "normal", bandwidth = 10)

# Load necessary library
library(MASS)  # For `isoMDS` function

# Generate example 3D coordinates
set.seed(123)
n <- 100
coordinates <- matrix(rnorm(n * 3), nrow = n, ncol = 3)

A_sp_nonspat <- read.csv("/Users/user/Documents/JHU/research/HRF_Est/Pain/A_sp_nonspat.csv")
A_sp_smooth <- A_sp_nonspat
numreg <- length(reglen)
cumsum_reglen <- c(0,cumsum(reglen))

for (i in 1:numreg){
  coordinates <- coord_list[[i]][1:reglen[i],]
  one_d_coordinates <- prcomp(coordinates, scale. = TRUE)$x[,1]
  smoothed_data <- ksmooth(one_d_coordinates, A_sp_nonspat[(cumsum_reglen[i]+1):cumsum_reglen[i+1],3], 
                           kernel = "normal", bandwidth = 0.2)$y
  if (length((cumsum_reglen[i]+1):cumsum_reglen[i+1]) >= 100){
    A_sp_smooth[(cumsum_reglen[i]+1):cumsum_reglen[i+1],3] <- smoothed_data
  }
  print(i)
}
write.csv(A_sp_smooth, file = "/Users/user/Documents/JHU/research/HRF_Est/Pain/A_sp_sm.csv", 
          row.names = FALSE)
 

# Residuals Analysis

library(forecast)
e_DD_sec1 <- e_DD_fit[1,1:300,]
e_DD_sec1_mean <- apply(e_DD_sec1,2,mean)
sec1 <- data.frame(t = 1:872, y = c(t(e_DD_sec1)))

ggplot(sec1) + 
  geom_line(aes(x=t, y = y, color = "Data"), alpha=0.2) +
  geom_line(data = data.frame(t = 1:872, y = e_DD_sec1_mean),
            aes(x = t, y = y, color = "Mean")) +
  scale_color_manual(name = "", values = c("Data" = "blue", "Mean" = "black"))+
  ggtitle("Residuals for Region 1: DD") +
  theme_bw()

arima_fit <- matrix(nrow = numreg, ncol = 3)
ind <- 1
for (i in 1:numreg){
  e_sp_sec1 <- apply(e_sp_fit[i,1:reglen[i],],2,mean)
  arima_fit[ind,] <- arimaorder(auto.arima(e_sp_sec1))
  ind <- ind + 1
  print(i)
}
par(mfrow=c(1,3))
hist(arima_fit[,1], xlab = "p", main="")
hist(arima_fit[,2], xlab = "q", main="")
hist(arima_fit[,3], xlab = "d", main="")


ts_e_DD_sec1 <- tsclean(ts(c(t(e_DD_sec1))))
pacf(ts_e_DD_sec1)

e_sp_sec1 <- e_sp_fit[1,1:300,]
e_sp_sec1_mean <- apply(e_sp_sec1,2,mean)
sec1 <- data.frame(t = 1:872, y = c(t(e_sp_sec1)))

ggplot(sec1) + 
  geom_line(aes(x=t, y = y, color = "Data"), alpha=0.2) +
  geom_line(data = data.frame(t = 1:872, y = e_sp_sec1_mean),
            aes(x = t, y = y, color = "Mean")) +
  scale_color_manual(name = "", values = c("Data" = "blue", "Mean" = "black"))+
  ggtitle("Residuals for Region 1: Spline") +
  theme_bw()


ts_e_sp_sec1 <- tsclean(ts(c(t(e_sp_sec1))))
pacf(ts_e_sp_sec1)

(fit <- Arima(ts_e_sp_sec1, order=c(4,1,0)))
checkresiduals(fit)


### Residuals average
e_DD_mean <- apply(e_DD_fit, c(1,2),mean)
e_DD_nonspat <- matrix(0,nrow = sum(reglen), ncol = 4)
ind<-1
for (i in 1:numreg){
  for (j in 1:reglen[i]){
    e_DD_nonspat[ind,1] <- e_DD_mean[i,j]
    e_DD_nonspat[ind,2:4] <- coord_list[[i]][j,]
    ind <- ind + 1
  }
}
