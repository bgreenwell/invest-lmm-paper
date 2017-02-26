########################################
#
# R code for manuscript titled "Linear Calibration with Random Coefficients via 
# the R Packages investr and lme4".
#
#  Author: Brandon M. Greenwell
# Journal: Journal of Statistical Software (JSS)
#    Year: 2014
#
########################################

# devtools::install_github("investr", username = "w108bmg")

# Preliminary -----------------------------------------------------------------

# Load packages
library(lattice)
library(nlme)
library(lme4)
library(rjags)
library(plyr)
library(investr)

# Load (original) bladder data 
# Bladder volume data
subject <- as.factor(rep(1:23, times = 8))
volume <- rep(c(10, 25, 50, 75, 100, 125, 150, 175), each = 23) / 10
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
bladder <- na.omit(data.frame(subject = subject, HD = HD, volume = volume))

# Scatterplot
xyplot(HD ~ volume, groups = subject, data = bladder, type = "b", alpha = 0.5)


# Fitted models ---------------------------------------------------------------

# Simple linear regression
lmFit <- lm(HD ~ volume, data = bladder)

# Random intercept and slope model
nlmeFit <- lme(HD ~ volume, random = list(subject = pdDiag(~volume)), 
               data = bladder)  # using nlme package
lme4Fit <- lmer(HD ~ volume + (0+1|subject) + (0+volume|subject), 
                data = bladder)  # using lme4 package

# Using invest function. Note that the CI for x0 based on the LMM is slightly
# shorter. This is not surprising since Var(Y0) = 17572.35 at x0.est for the 
# LMM, whereas, Var(Y0) is a constant 20249.88 in the LM.
x0.est <- invest(nlmeFit, y0 = 500, interval = "none")
invest(lmFit, y0 = 500)
invest(nlmeFit, y0 = 500)
invest(lmFit, y0 = 500, interval = "Wald")
invest(nlmeFit, y0 = 500, interval = "Wald")

# Simulation functions --------------------------------------------------------

# Calculate response variance at a particular predictor value
varY <- function(object, x) {
  vc <- as.data.frame(VarCorr(object))$sdcor
  vc[1]^2 + vc[2]^2*x^2 + vc[3]^2
}
sd.y0 <- sqrt(varY(lme4Fit, x = 8.015521))

# Generate data of different sample sizes
genData <- function(m, n) {
  volume <- rep(seq(from = 1, to = 17.5, length = n), times = m)
  subject <- rep(1:m, each = n)
  theta <- as.data.frame(VarCorr(lme4Fit))$sdcor
  beta <- getME(lme4Fit, "beta")
  alpha0 <- rnorm(m, mean = beta[1], sd = theta[1])
  alpha1 <- rnorm(m, mean = beta[2], sd = theta[2])
  HD <- rnorm(m*n, mean = alpha0[subject] + alpha1[subject]*volume,
              sd = theta[3])
  data.frame(subject, HD, volume)
}

# Fit model and calculate confidence intervals
getIntervals <- function(data, object, mean.response = FALSE) {
  eta <- if (mean.response) 500 else 500 + rnorm(1, mean = 0, sd = sd.y0)
  fit <- lme(HD ~ volume, random = list(subject = pdDiag(~volume)), 
             data = data, control = list(opt = "optim"))
  waldCI <- invest(fit, y0 = eta, interval = "Wald", tol = 1e-10,
                   mean.response = mean.response, lower = -10, upper = 30)
  invCI <- invest(fit, y0 = eta, interval = "inversion", tol = 1e-10,
                  mean.response = mean.response, lower = -10, upper = 30)
  c(lower.wald = waldCI$lower, upper.wald = waldCI$upper, 
    lower.inv = invCI$lower, upper.inv = invCI$upper)
}

# Function to summarize intervals
summarizeIntervals <- function(intervals) {
  res = cbind(length.wald = apply(intervals[, 1:2], 1, diff), 
              coverage.wald = ifelse(intervals[, 1] < 8.015521 & 
                                       intervals[, 2] > 8.015521, 1, 0), 
              length.inversion = apply(intervals[, 3:4], 1, diff), 
              coverage.inversion = ifelse(intervals[, 3] < 8.015521 & 
                                            intervals[, 4] > 8.015521, 1, 0))
  list(results = res, mean = apply(res, 2, mean), sd = apply(res, 2, sd))
}

# Simulation function
simFun <- function(x, nsim = 1000, seed = 101) {
  set.seed(seed)  # for reproducibility
  datasets <- rlply(nsim, genData(m = x[1], n = x[2]))
  cis <- ldply(datasets, getIntervals, object = nlmeFit, mean.response = FALSE)
  summarizeIntervals(cis)$mean
}

# Convert results to a single data frame
list2data <- function(x) {
  m <- n <- as.factor(c(5, 10, 30, 50, 100, 500))
  d <- expand.grid(n = n, m = m)
  z <- matrix(nrow = length(x), ncol = 4)
  for (i in seq(length(x))) z[i, ] <- x[[i]]
  z <- data.frame(z)
  z$m <- as.factor(d$m)
  z$n <- as.factor(d$n)
  setNames(z, c("length.wald", "coverage.wald", "length.inv", "coverage.inv",
                "m", "n"))
}

sample.sizes <- list("m = 5 and n = 5"   = c(5, 5), 
                     "m = 5 and n = 10"  = c(5, 10),
                     "m = 5 and n = 30"  = c(5, 30),
                     "m = 5 and n = 50"  = c(5, 50),
                     "m = 5 and n = 100" = c(5, 100),
                     "m = 5 and n = 500" = c(5, 500),
                     
                     "m = 10 and n = 5"   = c(10, 5), 
                     "m = 10 and n = 10"  = c(10, 10),
                     "m = 10 and n = 30"  = c(10, 30),
                     "m = 10 and n = 50"  = c(10, 50),
                     "m = 10 and n = 100" = c(10, 100),
                     "m = 10 and n = 500" = c(10, 500),
                     
                     "m = 30 and n = 5"   = c(30, 5), 
                     "m = 30 and n = 10"  = c(30, 10),
                     "m = 30 and n = 30"  = c(30, 30),
                     "m = 30 and n = 50"  = c(30, 50),
                     "m = 30 and n = 100" = c(30, 100),
                     "m = 30 and n = 500" = c(30, 500),
                     
                     "m = 50 and n = 5"   = c(50, 5), 
                     "m = 50 and n = 10"  = c(50, 10),
                     "m = 50 and n = 30"  = c(50, 30),
                     "m = 50 and n = 50"  = c(50, 50),
                     "m = 50 and n = 100" = c(50, 100),
                     "m = 50 and n = 500" = c(50, 500),
                     
                     "m = 100 and n = 5"   = c(100, 5), 
                     "m = 100 and n = 10"  = c(100, 10),
                     "m = 100 and n = 30"  = c(100, 30),
                     "m = 100 and n = 50"  = c(100, 50),
                     "m = 100 and n = 100" = c(100, 100),
                     "m = 100 and n = 500" = c(100, 500),
                     
                     "m = 500 and n = 5"   = c(500, 5), 
                     "m = 500 and n = 10"  = c(500, 10),
                     "m = 500 and n = 30"  = c(500, 30),
                     "m = 500 and n = 50"  = c(500, 50),
                     "m = 500 and n = 100" = c(500, 100),
                     "m = 500 and n = 500" = c(500, 500))
                     
# Sample run
library(parallel)
sim <- mclapply(sample.sizes, simFun, nsim = 1000, mc.cores = 6)
xyplot(coverage.inv ~ m|n, data = list2data(sim), type = "b", lwd = 2, pch = 19)

system.time(sim <- simFun(x = c(500, 500), nsim = 1000))

# Run simulations -------------------------------------------------------------
# library(doMC)
# registerDoMC(cores = 8)
# getDoParWorkers()  # check
sim.nlme.reg <- simulation(nsim = 1000, .progress = "text")
sim.nlme.cal <- simulation(nsim = 1000, mean.response = FALSE, .progress = "text")
save(sim.nlme.reg, sim.nlme.cal, 
     file = "/home/w108bmg/Desktop/Dropbox/sim_nlme.RData")

# library(doMC)
# registerDoMC(cores = 8)
nlmeSim <- simulation(nsim = 1000, mean.response = FALSE, .progress = "text")

# Plots for nlme model --------------------------------------------------------

# Plot coverage probabilities for nlme (regulation)
sim.nlme.cp <- list2data(sim.nlme)
levels(sim.nlme.cp$n) <- paste("n =", c(5, 10, 30, 50, 100, 500))
xyplot(cp ~ m|n, groups = method, data = sim.nlme.cp, type = "b", pch = 19, 
       xlab = "Number of subjects (m)", ylab = "Coverage probability",
       auto.key = list(columns = 2), 
       panel = function(x, y, ...) {
         panel.grid()
         panel.xyplot(x, y, ...)
         panel.abline(h = 0.95)
       })

# Plot coverage probabilities for nlme (regulation)
sim.nlme.cal.cp <- list2data(sim.nlme.cal)
levels(sim.nlme.cal.cp$n) <- paste("n =", c(5, 10, 30, 50, 100))
xyplot(cp ~ m|n, groups = method, data = sim.nlme.cal.cp, type = "b", pch = 19, 
       xlab = "Number of subjects (m)", ylab = "Coverage probability",
       auto.key = list(columns = 2), 
       panel = function(x, y, ...) {
         panel.grid()
         panel.xyplot(x, y, ...)
         panel.abline(h = 0.95)
       })


# # Plot coverage probabilities for lm (regulation)
# sim.lm.cp <- list2data(sim.lm)
# levels(sim.nlme.cp$n) <- paste("n =", c(5, 10, 30, 50, 100))
# xyplot(cp ~ m|n, groups = method, data = sim.lm.cp, type = "b", pch = 19, 
#        xlab = "Number of subjects (m)", ylab = "Coverage probability",
#        auto.key = list(columns = 2), 
#        panel = function(x, y, ...) {
#          panel.grid()
#          panel.xyplot(x, y, ...)
#          panel.abline(h = 0.95)
#        })
