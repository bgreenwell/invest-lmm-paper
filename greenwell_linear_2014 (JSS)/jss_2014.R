## Load packages ---------------------------------------------------------------
library(plyr)
library(lattice)
library(lme4)
library(nlme)
library(RColorBrewer); dark2 <- brewer.pal(8, "Dark2")
library(rjags)
library(investr)


## Load data -------------------------------------------------------------------
load("/home/w108bmg/Desktop/JSM 2014/bladder2.RData")

subject <- rep(1:23, times = 8)
volume <- rep(c(10, 25, 50, 75, 100, 125, 150, 175), each = 23)
HD <- c(13.2, 11.1, 10.3, NA, 4.8, 7.7, NA, 5.9, 1.9, 6.5, 19.8, 
        14.6, NA, NA, 9.7, 17.2, 10.6, 19.3, 8.5, 6.9, 8.1, 14.8, 13.7, 
        27.4, 27.5, 15, 10, 18.6, 12.6, 24, 28.4, 12.5, 16.7, 29.6, 
        27.1, 14, 18.7, 20.3, 35.8, 23.6, 37.4, 31.3, 23.7, 22, 34.3, 
        28.5, 41.6, 58.1, 34.2, 28.8, 29.9, 31.4, 46.9, 44.4, 26.8, 
        30.6, 51.7, 49.8, 19.1, 35.8, 38.9, 41.4, 49.9, 58.6, 54.8, 44, 
        39.1, 58.5, 41.5, 60.1, 78.8, 49.4, 46.4, 39.4, 45.3, 50.4, 
        70.7, 54.4, 41.8, 72.2, 67.5, 39.2, 49.6, 65.1, 69.7, 67.7, 
        73.7, 78.3, 65.7, 44.7, 72.1, 59.8, 73.9, 91.5, 71.3, 54.8, NA, 
        48, 67.8, 89.4, 63.1, 49.6, 81.9, 79.1, 48.7, 65.6, 65.1, 81.9,
        87.7, 79.4, 93, 80.3, 68.9, 90.9, 77.5, 85.5, 98.3, 81.3, 69.4, 
        NA, 66.6, 81, 105.8, 83.5, 60.8, 95.1, 95.1, 67, 85.3, 86.9, 
        96.6, 89.3, 102.6, NA, 93.6, 93.3, 105, 92.9, 95.6, 111.4, 94, 
        73.9, NA, NA, 91.2, 113.5, 114.5, 80.1, 115.4, 109.8, 72.7, 
        90.4, 98.6, 115, 108, 110.9, NA, 99.2, 102.4, 117.5, 99.4, 
        107.4, 121, 104.3, NA, NA, NA, 99.8, 127.3, 124, 87.1, NA, NA, 
        NA, NA, 107.2, 117, 114.8, 122.4, NA, 112.2, 104.7, 124.2, 113)
bladder <- data.frame(subject = subject, HD = HD, volume = volume)
bladder <- na.omit(bladder)
bladder$volume <- bladder$volume/10

## Figures ---------------------------------------------------------------------

## Figure path
fig.path <- "/home/w108bmg/Desktop/greenwell-linear-2014/"

## Figure 1
# algorithm box

## Figure 2
pdf(file = paste(fig.path, "bladder-spaghettiplots.pdf", sep = "/"), width = 7, height = 4)
p1 <- xyplot(HD ~ volume, groups = subject, data = bladder, type = "l",
             xlab = "Volume (cl)", ylab = "Height (mm) times Depth (mm)")
p2 <- xyplot(HD ~ volume, groups = subject, data = bladder2, type = "l",
             xlab = "Volume (cl)", ylab = "Height (mm) times Depth (mm)")
print(p1, position = c(0, 0, 0.5, 1), more=TRUE)
print(p2, position = c(0.5, 0, 1, 1))
dev.off()

## LMM for bladder volume data using nlme package
bladder.nlme <- lme(HD ~ volume + I(volume^2), data = bladder,
                    random = list(subject = pdDiag(~volume)))

## LMM for bladder volume data using lme4 package
bladder.lme4 <- lmer(HD ~ volume + I(volume^2) + (0+1|subject) + 
                       (0+volume|subject), data = bladder)

## LMM for transformed bladder volume data using nlme package
bladder2.nlme <- lme(HD ~ volume, random = list(subject = pdDiag(~volume)),
                     data = bladder2)

## LMM for transformed bladder volume data using lme4 package
bladder2.lme4 <- lmer(HD ~ volume + (0+1|subject) + (0+volume|subject),
                      data = bladder2)

## Generate new observations
var.y <- function(object, x) {
  vc <- as.data.frame(VarCorr(object))$sdcor
  vc[1]^2 + vc[2]^2*x^2 + vc[3]^2
}
sd.y0 <- sqrt(var.y(bladder2.lme4, x = 8.015521))



## Simulation study ------------------------------------------------------------

## Simulation data
genData <- function(nsim = 1) {
  bladder2.sim <- vector("list", length = nsim)
  for (i in 1:nsim) {
    bladder2.sim[[i]] <- bladder2
    bladder2.sim[[i]]$HD <- simulate(bladder2.lme4)[[1]]
    names(bladder2.sim[[i]]) <- names(bladder2)
  }
  bladder2.sim
}

## Function to generate fake bladder data with different sample sizes
bladder.fake <- function(m, n) {
  volume <- rep(seq(from = 1, to = 17.5, length = n), times = m)
  subject <- rep(1:m, each = n)
  theta <- as.data.frame(VarCorr(bladder.lme4))$sdcor
  beta <- getME(bladder.lme4, "beta")
  alpha0 <- rnorm(m, mean = beta[1], sd = theta[1])
  alpha1 <- rnorm(m, mean = beta[2], sd = theta[2])
  HD <- rnorm(m*n, mean = alpha0[subject] + alpha1[subject]*volume + 
                beta[3]*volume^2, sd = theta[3])
  data.frame(subject, HD, volume)
}

## Function to generate fake transformed bladder data with different sample 
## sizes
bladder2.fake <- function(m, n) {
  volume <- rep(seq(from = 1, to = 17.5, length = n), times = m)
  subject <- rep(1:m, each = n)
  theta <- as.data.frame(VarCorr(bladder2.lme4))$sdcor
  beta <- getME(bladder2.lme4, "beta")
  alpha0 <- rnorm(m, mean = beta[1], sd = theta[1])
  alpha1 <- rnorm(m, mean = beta[2], sd = theta[2])
  HD <- rnorm(m*n, mean = alpha0[subject] + alpha1[subject]*volume, 
              sd = theta[3])
  data.frame(subject, HD, volume)
}

## Function to generate a list of simulated data frames
makeDataList <- function(nsim, m, n, linear = TRUE) {
  if (linear) {
    rlply(nsim, bladder2.fake(m, n)) # simulate from transformed data
  } else {
    rlply(nsim, bladder.fake(m, n)) # simulate from original data
  }
}

## Function to fit model and calculate classical intervals for data
getIntervals <- function(d) {
  fm <- update(bladder2.nlme, data = d, control = list(opt = "optim"))
#                control = list(maxItr = 500, msMaxIter = 500, msMaxEval = 500,
#                               niterEM = 100))
  res1 <- invest(fm, y0 = 500, interval = "Wald", mean.response = T)
  res2 <- invest(fm, y0 = 500, interval = "inversion", mean.response = T,
                 lower = 0, upper = 25)
  c(lower.W = res1$lower, upper.W = res1$upper, 
    lower.I = res2$lower, upper.I = res2$upper)
#   fm <- update(bladder2.nlme, data = d, control = lmeControl(opt = "optim"))
#   Y0 <- 500 + rnorm(1, mean = 0, sd = sd.y0) # FIXME: Is this the best way?
#   res1 <- invest(fm, y0 = Y0, interval = "Wald")
#   res2 <- invest(fm, y0 = Y0, interval = "inversion", lower = -10, upper = 40)
#   c(lower.W = res1$lower, upper.W = res1$upper, 
#     lower.I = res2$lower, upper.I = res2$upper)
}

## Function to summarize intervals
summarizeIntervals <- function(intervals) {
  length.W <- apply(intervals[, 1:2], 1, diff)
  length.I <- apply(intervals[, 3:4], 1, diff)
  cover.W <- ifelse(intervals[, 1] < 8.015521 & intervals[, 2] > 8.015521, 1, 0)
  cover.I <- ifelse(intervals[, 3] < 8.015521 & intervals[, 4] > 8.015521, 1, 0)
  res = cbind(length.W = length.W, coverage.W = cover.W, 
              length.I = length.I, coverage.I = cover.I)
  list(results = res, summary = rbind(mean = apply(res, 2, mean),
                                        sd = apply(res, 2, sd)))
}

## Row one of simulation results table -----------------------------------------
set.seed(111)
z <- makeDataList(1000, m = 5, n = 5)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(112)
z <- makeDataList(1000, m = 5, n = 10)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(113)
z <- makeDataList(1000, m = 5, n = 30)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(114)
z <- makeDataList(1000, m = 5, n = 50)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(115)
z <- makeDataList(1000, m = 5, n = 50)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

## Row two of simulation results table -----------------------------------------
set.seed(121)
z <- makeDataList(1000, m = 10, n = 5)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(122)
z <- makeDataList(1000, m = 10, n = 10)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(123)
z <- makeDataList(1000, m = 10, n = 30)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(124)
z <- makeDataList(1000, m = 10, n = 50)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(125)
z <- makeDataList(1000, m = 10, n = 100)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

## Row three of simulation results table ---------------------------------------
set.seed(131)
z <- makeDataList(1000, m = 30, n = 5)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(132)
z <- makeDataList(1000, m = 30, n = 10)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(133)
z <- makeDataList(1000, m = 30, n = 30)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(134)
z <- makeDataList(1000, m = 30, n = 50)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(135)
z <- makeDataList(1000, m = 30, n = 100)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

## Row four of simulation results table ----------------------------------------
set.seed(141)
z <- makeDataList(1000, m = 50, n = 5)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(142)
z <- makeDataList(1000, m = 50, n = 10)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(143)
z <- makeDataList(1000, m = 50, n = 30)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(144)
z <- makeDataList(1000, m = 50, n = 50)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(145)
z <- makeDataList(1000, m = 50, n = 100)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

## Row five of simulation results table ----------------------------------------
set.seed(151)
z <- makeDataList(1000, m = 100, n = 5)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(152)
z <- makeDataList(1000, m = 100, n = 10)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(153)
z <- makeDataList(1000, m = 100, n = 30)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(154)
z <- makeDataList(1000, m = 100, n = 50)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)

set.seed(155)
z <- makeDataList(1000, m = 100, n = 100)
intervals <- ldply(z, getIntervals, .progress = "text")
round(summarizeIntervals(intervals)$summary, 2)
