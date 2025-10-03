#### 1. Import libraries ####
library(reshape2)
library(gdata)
library(DescTools)
library(R.matlab)
library(data.table)
library(zoo)
library(MTS)

#### 2. Paths ####
out_path <- ""
ts_input <- ""
regions <- ""
hcp.dir <- ""



## Background

# VAR analysis of fMRI data on healthy participants. Very briefly, p = 86,
# T = 1200, TR= 720 ms for each subject.  We focus on p = 6,
# T = 1200 for each subject. 


## Data cleaning


#### 3. Get subject id #### 
fc_subjects <- sub("^(.*?)_.*$", "\\1", 
                   list.files(path = ts_input)) 
sc_subjects <- as.character(dat.SC$subject)
subjects <- intersect(fc_subjects, sc_subjects)
length(subjects)

#### 4. Extract FC and SC #### 
n = length(subjects)
r = length(regions) 
tp = 1200
TS <- array(0, c(n,tp,r))
SC <- array(0, c(n,r,r))

for (i in 1:n) {
  id <- subjects[i]
  # print(id)
  
  # Generate FC
  load(paste0(ts_input, id, "_TS.RData"))
  # 4 replicates 4 x 1200 time points = 4800
  ts_list <- list()
  for (j in 1:4) {
    reps <- nrow(timeseries)/4
    ts <- timeseries[((j-1)*reps + 1):(j*reps), 1:r] 
    if (all(is.na(ts))) {
      next
    } else {
      # pearson cor
      ts_list[[j]] <- ts
    }
  }
  
  # remove null replicate
  ts_list[sapply(ts_list, is.null)] <- NULL
  
  
  if (length(ts_list) == 4) {
    # average 4 ts 
    TS[i,,] <- as.matrix(ts_list[[1]] + ts_list[[2]] + 
                           ts_list[[3]] +ts_list[[4]]) / 4 
  } else {
    next
  }
  
  # Get SC sift2volnorm version
  SC[i,,] <- as.array(dat.SC$sift2volnorm[which(dat.SC$subject==as.numeric(id)), ,])
}

dim(TS)

# Combine all regions
regions <- data.table(r_num = 1:r, r_name = regions)


## Focus on 6 subcortical regions

subcort <- c("Thalamus-Proper", "Hippocampus", "Amygdala")
subcort_indx <- which(grepl(paste(subcort, collapse="|"), regions$r_name))
subcort_r <- regions[subcort_indx, r_name]
TS_subcort <- TS[, , subcort_indx]
SC_subcort <- SC[, , subcort_indx]


## Plot TS of 1 subject

# All times series

i = 1
x <- data.frame(TS_subcort[i, ,])
dim(x)
names(x) <- c("L-Thalamus", "L-Hippo", "L-Amyg",
              "R-Thalamus", "R-Hippo", "R-Amyg")
head(x)

# normalize the data
x <- scale(x)

par(mfrow = c(1,2))
matplot(x, type = 'l', xlab="Time", ylab="")
# ts.plot(x)
plot.ts(x, main="")
par(mfrow = c(1,1))



# ## ACF and PACF
# L- and R-Thalamus: ACF shows exponential decay, PACF shows significant up to lag 7
# L-Hippo: ACF also shows exponential decay, slightly significant at the end, PACF shows significant at 1, 2, 27, 29
# R-Hippo: ACF shows exponential decay, PACF shows significant up to lag 3
# L-Amyg: ACF shows exponential decay, PACF shows significant up to lag 7
# R-Amyg: ACF shows exponential decay, PACF shows significant up to lag 3


for (j in 1:6) {
  par(mfrow = c(1,2))
  acf(x[, j])
  pacf(x[, j])
  par(mfrow = c(1,1))
}

## Cross correlation

# Separate into each hemisphere


left_x <- x[, 1:3]

matplot(left_x, type = 'l', xlab="Time", ylab="")
ts.plot(left_x)
names(left_x) <- gsub("Left-", "", names(left_x))
plot.ts(left_x, main="")

acf(left_x)              # cross-corrleation function 2x2
acf(left_x[, 1:2])              # cross-corrleation function 2x2
acf(left_x[, 2:3])              # cross-corrleation function 2x2
acf(left_x[, c(1,3)])              # cross-corrleation function 2x2


right_x <- x[, 4:6]

matplot(right_x, type = 'l', xlab="Time", ylab="")
ts.plot(right_x)
names(right_x) <- gsub("Right-", "", names(right_x))
plot.ts(right_x, main="")

acf(right_x)              # cross-corrleation function 2x2
acf(right_x[, 1:2])              # cross-corrleation function 2x2
acf(right_x[, 2:3])              # cross-corrleation function 2x2
acf(right_x[, c(1,3)])              # cross-corrleation function 2x2


# Same region in each hemisphere

thal_x <- x[, c(1,4)]
acf(thal_x)              # cross-corrleation function 2x2

hippo_x <- x[, c(2,5)]
acf(hippo_x)              # cross-corrleation function 2x2

amyg_x <- x[, c(3,6)]
acf(amyg_x)              # cross-corrleation function 2x2


# ## VAR model
# Ljung-Box test 
# Note that it is applied to the residuals of a fitted ARIMA model, not the original series, and in such applications the hypothesis actually being tested is that the residuals from the ARIMA model have no autocorrelation.
# Use Augmented Dickey-Fuller Test

# mq(x, lag = 5)    
library(fpp)
for (j in 1:6) {
  print(adf.test(x[,j]))
}


# VAR(p) model

y <- as.ts(x)  

#acf(y)

fit.VARp <- VARorder(y, maxp = 10, output = T)

plot(1:11, fit.VARp$aic, type = 'l')
points(which.min(fit.VARp$aic)-1, fit.VARp$aic[which.min(fit.VARp$aic)], pch = 20, col = 2, cex = 3)
which.min(fit.VARp$aic) # (subtract 1 for lag order since lag 0 is included) 



# Fit a VAR(2) model 
VAR(y, p = 2, output = T, include.mean = F, fixed = NULL)
fit.VAR2 = VAR(y, p = 2, output = T, include.mean = F, fixed = NULL)
names(fit.VAR2)

fit.VAR2$coef
fit.VAR2$Phi
fit.VAR2$Ph0
fit.VAR2$Sigma

# se's for regression parameters:
fit.VAR2$secoef

# t-ratios for regression parameters:
t = fit.VAR2$coef / fit.VAR2$secoef

# Estimated residual series
a.hat <- fit.VAR2$residuals  
head(a.hat)
nrow(y)
nrow(a.hat) # removed first 2 rows since we fit VAR(2)


# Diagnostic checking
par(mfrow = c(2,2))
ccf(a.hat[,1],a.hat[,2])
matplot(a.hat, type = 'l')
plot(a.hat[,1], type = 'l') ; abline(h = 0)
plot(a.hat[,2], type = 'l') ; abline(h = 0)
par(mfrow = c(1,1))


plot(a.hat[,1],a.hat[,2])
plot(a.hat[,1],a.hat[,4])


par(mfrow = c(2,2))
hist(a.hat[,1])
hist(a.hat[,2])
hist(a.hat[,3])

hist(a.hat[,4])
hist(a.hat[,5])
hist(a.hat[,6])


qqnorm(a.hat[,1]) ; qqline(a.hat[,1])
qqnorm(a.hat[,2]) ; qqline(a.hat[,2])
qqnorm(a.hat[,3]) ; qqline(a.hat[,3])

qqnorm(a.hat[,4]) ; qqline(a.hat[,4])
qqnorm(a.hat[,5]) ; qqline(a.hat[,5])
qqnorm(a.hat[,6]) ; qqline(a.hat[,6])

par(mfrow = c(1,1))


acf(a.hat)

mq(a.hat, lag = 5)    

grangertest(a.hat[,1], a.hat[,4], order = 2)
grangertest(a.hat[,4], a.hat[,1], order = 2)

grangertest(a.hat[,2], a.hat[,5], order = 2)
grangertest(a.hat[,5], a.hat[,2], order = 2)

grangertest(a.hat[,3], a.hat[,6], order = 2)
grangertest(a.hat[,6], a.hat[,3], order = 2)




# Multivariate Ljung-Box Q(m) statistics w and w/o adjustment
# g <- p * k^2 
g <- 2*2^2 ; g  # number of regression parameters (ignoring intercept)
mq(a.hat,10)
mq(a.hat,10, adj = g)
mq(a.hat,15)
mq(a.hat,15, adj = g)

# MTSdiag(fit.VAR2, gof = 24, adj = 0, level = F)

# acf(a.hat^2)
MarchTest(a.hat, lag = 10)


# Refine a fitted VAR model by removing simultaneously insignificant parameters
refVAR(fit.VAR2, fixed = NULL, thres = 2)


## Fit model for all participants 


ar1 <- array(0, c(n,6,6))
ar2 <- array(0, c(n,6,6))

for (i in 1:49) {
  x <- data.frame(TS_subcort[i, ,])
  names(x) <- c("L-Thalamus", "L-Hippo", "L-Amyg",
                "R-Thalamus", "R-Hippo", "R-Amyg")
  
  # normalize the data
  x <- scale(x)
  
  y <- as.ts(x)
  fit.VAR2 = MTS::VAR(y, p = 2, include.mean = F, fixed = NULL, output=F)
  refit.VAR2 = MTS::refVAR(fit.VAR2, fixed = NULL, thres = 2)
  ar1[i, , ] <- t(refit.VAR2$coef[1:6,])
  ar2[i, , ] <- t(refit.VAR2$coef[7:12,])
  
}

# Correlation of ar1 with SC


cor(scale((SC_subcort[, 4,1])), ar1[, 4, 1])
cor(scale((SC_subcort[, 5,2])), ar1[, 5, 2])
cor(scale((SC_subcort[, 6,3])), ar1[, 6, 3])

par(mfrow = c(3,1))
plot(scale((SC_subcort[, 4,1])), ar1[, 4, 1],
     ylab="Thalamus",
     xlab="", xlim=c(-3,3), ylim=c(0,0.5))
plot(scale((SC_subcort[, 5,2])), ar1[, 5, 2],
     ylab="Hippocampus",
     xlab="", xlim=c(-3,3), ylim=c(0,0.5))
plot(scale((SC_subcort[, 6,3])), ar1[, 6, 3],
     ylab="Amygdala",
     xlab="Scaled structural connectivity",
     xlim=c(-3,3), ylim=c(0,0.5))
par(mfrow = c(1,1))



# Reverse


cor(scale(log(SC_subcort[, 4,1])), ar1[, 1, 4])
cor(scale(log(SC_subcort[, 5,2])), ar1[, 2, 5])
cor(scale(log(SC_subcort[, 6,3])), ar1[, 3, 6])

par(mfrow = c(3,1))
plot(scale((SC_subcort[, 4,1])), ar1[, 1, 4],
     ylab="Thalamus",
     xlab="", xlim=c(-3,3), ylim=c(0,0.5))
plot(scale((SC_subcort[, 5,2])), ar1[, 2, 5],
     ylab="Hippocampus",
     xlab="", xlim=c(-3,3), ylim=c(0,0.5))
plot(scale((SC_subcort[, 6,3])), ar1[, 3, 6],
     ylab="Amygdala",
     xlab="Scaled log structural connectivity",
     xlim=c(-3,3), ylim=c(0,0.5))
par(mfrow = c(1,1))



# Correlation of ar2 with SC


cor(scale((SC_subcort[, 4,1])), ar2[, 4, 1])
cor(scale((SC_subcort[, 5,2])), ar2[, 5, 2])
cor(scale((SC_subcort[, 6,3])), ar2[, 6, 3])

par(mfrow = c(3,1))
plot(scale((SC_subcort[, 4,1])), ar2[, 4, 1],
     ylab="Thalamus",
     xlab="", xlim=c(-3,3), ylim=c(0,0.5))
plot(scale((SC_subcort[, 5,2])), ar2[, 5, 2],
     ylab="Hippocampus",
     xlab="", xlim=c(-3,3), ylim=c(0,0.5))
plot(scale((SC_subcort[, 6,3])), ar2[, 6, 3],
     ylab="Amygdala",
     xlab="Scaled structural connectivity",
     xlim=c(-3,3), ylim=c(0,0.5))
par(mfrow = c(1,1))

