#Preprocessing Space Top Ex Data

# library
library(readr)
library(tidyverse)
library(scatterplot3d)
library(ggplot2)
library(stats)
library(pracma)
library(scatterplot3d)
library(graphics)
setwd("/Users/user/Dropbox (Personal)/SpaceTopEx/subject-0002")


TR <- 0.460
nvol <- Tlen <- 872

cue_1 <- round(read_csv("EVs/Cue-Run1.csv")[,-3]/TR)
exp_1 <- round(read_csv("EVs/Expectrating-Run1.csv")[,-3]/TR)
stim_1 <- read_csv("EVs/Stimulus-Run1.csv")[,-(3:4)]
out_1 <- round(read_csv("EVs/Outcomerating-Run1.csv")[,-3]/TR)

stim_1$duration <- stim_1$duration - stim_1$onset
stim_1 <- round(stim_1/TR)

exp_1$duration[is.na(exp_1$duration)] <- round(mean(exp_1$duration[!is.na(exp_1$duration)]))

#cue
S1 <- S2 <- S3 <- S4 <- rep(0, nvol)
for (i in 1:nrow(cue_1)) {
  S1[cue_1[i,1]:(cue_1[i,1] + cue_1[i,2])] <- 1
}

#expect rating
for (i in 1:nrow(exp_1)) {
  S2[exp_1[i,1]:(exp_1[i,1] + exp_1[i,2])] <- 1
}

#stimulus
for (i in 1:nrow(stim_1)) {
  S3[stim_1[i,1]:(stim_1[i,1] + stim_1[i,2])] <- 1
}

#outcome rating
for (i in 1:nrow(out_1)) {
  S4[out_1[i,1]:(out_1[i,1] + out_1[i,2])] <- 1
}

Runc <- list(S1, S2, S3, S4)
save(Runc, file = "Runc.RData")

# Plot stimuli function for SpaceTopEx

df_new <- data.frame(
  x = rep(1:length(Runc[[1]]), length(Runc)),
  y = unlist(Runc),
  S = factor(rep(c("cue", "expect", "stimulus", "outcome"), each = length(Runc[[1]])))
)

ggplot(df_new, aes(x = x, y = y, group = S)) +
  geom_line(aes(color = S)) +
  labs(title = "Stimuli function ", x = "t", y = "") 

# Remove nuisance variables
X <- as.matrix(read.table("EVs/Nuisance-run1.txt", sep = ",",
                 header = FALSE, quote = "", stringsAsFactors = FALSE))


folders <- list.files("/Users/user/Dropbox (Personal)/SpaceTopEx/subject-0002/RegionData", 
                     recursive = FALSE)

numreg <- 268
dat_list <- vector("list", length = numreg)
coord_list <- vector("list", length = numreg)

setwd("/Users/user/Dropbox (Personal)/SpaceTopEx/subject-0002/RegionData")
for (k in 1:numreg){
  datastr <- folders[2*k]
  coordstr <- folders[2*k-1]
  
  # fmri timecourse data
  raw_t1 <- t(read.table(datastr, sep = ",",
                         header = FALSE, quote = "", stringsAsFactors = FALSE))
  coord_t1 <- t(read.table(coordstr, sep = ",",
                         header = FALSE, quote = "", stringsAsFactors = FALSE))
  dat1 <- matrix(0,nrow(raw_t1), ncol(raw_t1))
  
  # timecourse error term (- projection on regressors)
  for (i in 1:nrow(raw_t1)) {
    y <- matrix(raw_t1[i,], ncol = 1)
    e <- y - X %*% base::solve(t(X) %*% X, t(X)) %*% y
    dat1[i,] <- t(e)
  }
  dat_list[[k]] <- dat1
  coord_list[[k]] <- coord_t1
}
save(dat_list, coord_list, file = "data_list.RData")



num <- nrow(dat1)

# Plot average timecourse
df_tc <- data.frame(t = rep(1:nvol, 300), y = c(t(dat1)), voxel = rep(1:300, each = nvol))

df_tc_tmp <- df_tc 
ggplot(data = df_tc, aes(x = t, group = voxel)) + 
  geom_line(aes(y = y)) + 
  xlab("t") + 
  theme_bw() +
  ggtitle("Distribution of Timecourse Data For Region 1")

