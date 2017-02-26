## Simulation study for JStatSoft paper

## Preliminary -----------------------------------------------------------------

## Load packages
library(investr)  # inverse estimation functions
library(nlme)     # LMM functions (S3)
library(lme4)     # LMM functions (S4 and parametric bootstrap)
library(plyr)     # utility functions
library(lattice)

## Bladder data
subject <- rep(1:23, times = 8)
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
bladder <- data.frame(subject = subject, HD = HD, volume = volume)
bladder <- na.omit(bladder)

xyplot(HD^(3/2) ~ volume, groups = subject, data = bladder, type = "b")

## Random intercept and slope model (lme4 version)
bladder_lme4 <- lmer(HD^(3/2) ~ volume + (0+1|subject) + (0+volume|subject), 
                     data = bladder)

bladder_nlme <- lme(HD^(3/2) ~ volume, random = list(subject = pdDiag(~volume)), 
                    data = bladder) 

varY <- function(object, x) {
  vc <- as.data.frame(VarCorr(object))$sdcor
  vc[1]^2 + vc[2]^2*x^2 + vc[3]^2
}
sd_y0 <- sqrt(varY(bladder_lme4, x = 8.015521))

## Functions -------------------------------------------------------------------

## Simulate transformed bladder data
genData <- function(m = 23, n = 8) {
  volume <- rep(seq(from = 1, to = 17.5, length = n), times = m)
  subject <- rep(1:m, each = n)
  theta <- as.data.frame(VarCorr(bladder_lme4))$sdcor
  beta <- getME(bladder_lme4, "beta")
  alpha0 <- rnorm(m, mean = beta[1], sd = theta[1])
  alpha1 <- rnorm(m, mean = beta[2], sd = theta[2])
  HD <- rnorm(m*n, mean = alpha0[subject] + alpha1[subject]*volume, 
              sd = theta[3])
  data.frame(subject, HD, volume)
}

## Replicate simulation
simulate <- function(nsim = 10, m = 23, n = 8, y0 = 500, sd.y0 = 132.5607,
                     mean.response = FALSE, ...) {
  rdply(nsim, {
    Y0 <- if (mean.response) y0 else y0 + rnorm(1, mean = 0, sd = sd.y0)
    fm <- lme(HD ~ volume, random = list(subject = pdDiag(~volume)), 
              data = genData(m = m, n = n),
              control = lmeControl(opt = "optim"))
    w <- invest(fm, y0 = Y0, interval = "Wald", mean.response = mean.response, 
                tol= 1e-10)
    i <- invest(fm, y0 = Y0, interval = "inversion", 
                mean.response = mean.response, lower = -5, upper = 40, 
                tol= 1e-10)
    c("w_est" = w$estimate, "w_lwr" = w$lower, "w_upr" = w$upper,
      "i_est" = i$estimate, "i_lwr" = i$lower, "i_upr" = i$upper)
  }, ...)
}

## Summarize simmulation
simSummary <- function(x, x0 = 8.015521) {
  c("w_cp"  = mean(ifelse(x$w_lwr < x0 & x$w_upr > x0, 1, 0)),
    "i_cp"  = mean(ifelse(x$i_lwr < x0 & x$i_upr > x0, 1, 0)),
    "w_len" = mean(x$w_upr - x$w_lwr),
    "i_len" = mean(x$i_upr - x$i_lwr))
} 

## Data frame of sample sizes
m <- n <- c(5, 10, 15, 20, 25, 30, 50, 100)
grid <- setNames(expand.grid(m, n), c("m", "n"))

## Run simulation (calibration)
system.time({
set.seed(101)
res <- vector("list", length = nrow(grid))
for (i in seq_len(nrow(grid))) {
  res[[i]] <- simulate(nsim = 1000, m = grid[i, 1], n = grid[i, 2])
  print(paste(round(i/length(res) * 100), "% complete", sep = ""))
}
})
df_wide <- cbind(grid, ldply(res, simSummary))

## Plot results
library(reshape2)
df_coverage <- melt(z, id.vars = c("m", "n"), measure.vars = c("w_cp", "i_cp"),
                    variable.name = "interval", value.name = "coverage")
xtabs(coverage ~ m + n + interval, data = df_coverage)
xyplot(coverage ~ as.factor(m) | as.factor(n), groups = interval, type = "b", 
       data = df_coverage,
       panel = function(...) {
         panel.abline(h = 0.95, col = "black")
         panel.xyplot(...)
})

# 45.417 