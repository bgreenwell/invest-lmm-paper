## Load packages ###############################################################
library(lattice) # for plotting
library(car)     # for deltaMethod() function
library(nlme)    # for fitting LMMs
library(lme4)    # for fitting LMMs
library(boot)    # for calculating bootstrap CI's


## Set up paths for figures ####################################################
work.dir <- "I://setup//Desktop//Linear calibration with grouped data - Biometrics"
home.dir <- "/home/w108bmg/Desktop/Linear calibration with grouped data - Biometrics"
desk.dir <- NULL


## Bladder data ################################################################
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
Bladder <- groupedData(HD ~ volume | subject, data = bladder)

## Fit models ##################################################################
bladder.nlme <- lme(HD ~ volume + I(volume^2), data = Bladder,
                    random = list(subject = pdDiag(~volume))) 
         
         
## Point estimate ##############################################################
y0 <- 70
x0.est <- uniroot(function(x) {
  predict(bladder.nlme, list(volume = x), level = 0) - y0 }, 
  interval = range(bladder$volume), tol = 1e-10, 
  maxiter = 1000)$root


## Wald-based interval #########################################################
b <- as.numeric(fixef(bladder.nlme))
var.y0 <- getVarCov(bladder.nlme)[1, 1] +
  x0.est^2*getVarCov(bladder.nlme)[2, 2] + 
  summary(bladder.nlme)$sigma^2
covmat <- diag(4)
covmat[1:3, 1:3] <- vcov(bladder.nlme) 
covmat[4, 4] <- var.y0
params <- c(b0 = b[1], b1 = b[2], b2 = b[3], Y0 = y0)
g <- "(-b1 + sqrt(b1^2 - 4*b2*(b0-Y0))) / (2*b2)"
dm <- deltaMethod(params, g = g, vcov. = covmat)
wald.ci <- c(x0.est, x0.est + qnorm(c(0.025, 0.975))*dm$SE)


## Inversion interval ##########################################################
predFun <- function(x) {
  z <- list(x)
  names(z) <- "volume"
  fit <- predict(bladder.nlme, newdata = z, level = 0)
  se.fit <- sqrt(diag(cbind(1, unlist(z), unlist(z)^2) %*%
                        bladder.nlme$varFix %*%
                        t(cbind(1, unlist(z), unlist(z)^2))))
  list(fit = fit, se.fit = se.fit)
}
bounds <- function(x) {
  z <- list(x)
  names(z) <- "volume"
  pred <- predFun(x)
  (y0 - pred$fit)^2/(var.y0 + pred$se.fit^2) - qnorm(0.975)^2
}
lower <- uniroot(bounds, interval = c(min(bladder$volume), x0.est),
    tol = 1e-10, maxiter = 1000)$root
upper <- uniroot(bounds, interval = c(x0.est, max(bladder$volume)),
    tol = 1e-10, maxiter = 1000)$root
inversion.ci <- c(x0.est, c(lower, upper))

## Compare with invest function
inversion.ci
invest(bladder.nlme, y0 = 70, interval = "inversion")
wald.ci
invest(bladder.nlme, y0 = 70, interval = "Wald")

## Parametric bootstrap ########################################################

## Refit model using lme4 package
bladder.lme4 <- lmer(HD ~ volume + I(volume^2) + (0+1|subject) +
  (0+volume|subject), data = bladder)

## Calculate standard deviation of Y0
sd.y0 <- sqrt(VarCorr(bladder.lme4)[[1]][1] +
  x0.est^2*VarCorr(bladder.lme4)[[2]][1] + sigma(bladder.lme4)^2)

## Bootstrap function
bootFun <- function(.) {
  ## Bootstrap point estimate directly
  beta.boot <- as.numeric(fixef(.)) # fixed effects
  if (all(getME(., "y") == bladder$HD)) {
    y0.boot <- y0 # return original estimate
  } else {
    y0.boot <- y0 + rnorm(1, mean = 0, sd = sd.y0)
  }
  x0.boot <- (-beta.boot[2] + sqrt(beta.boot[2]^2 - 
    4*beta.boot[3]*(beta.boot[1]-y0.boot)))/(2*beta.boot[3])
                                   
  ## Adjusted inversion interval
  covb <- vcov(.)               # var-cov matrix for beta
  v1.boot <- VarCorr(.)[[1]][1] # variance for random intercept 
  v2.boot <- VarCorr(.)[[2]][1] # variance for random slope
  v3.boot <- sigma(.)^2         # residual variance
  mu0.boot <- as.numeric(crossprod(beta.boot, 
    c(1, x0.est, x0.est^2)))
  var.y0 <- v1.boot + x0.est^2*v2.boot + v3.boot
  var.mu0 <- t(c(1, x0.est, x0.est^2)) %*% covb %*% 
    c(1, x0.est, x0.est^2)
  Q.boot <- as.numeric((y0.boot - mu0.boot) / 
    sqrt(var.y0 + var.mu0))
  
  ## Return results
  c(x0.boot, Q.boot)
}

## Run bootstrap simulation
# set.seed(0101)
# (bladder.pb <- bootMer(bladder.lme4, FUN = bootFun, nsim = 9999))
# save(bladder.pb, file = paste(home.dir, "bladderPB.RData", sep="/"))
load(paste(home.dir, "bladderPB.RData", sep="/"))

## Calculate bootstrap intervals for x0
boot.ci(bladder.pb, index = 1, type = c("norm", "perc"))

## Bootstrap adjusted inversion interval
Q.025 <- quantile(as.numeric(bladder.pb$t[, 2]), 0.025)
Q.975 <- quantile(as.numeric(bladder.pb$t[, 2]), 0.975)
upper <- function(x) {
  z <- list(x)
  names(z) <- "volume"
  pred <- predFun(x)
  (y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - Q.975
}
lower <- function(x) {
  z <- list(x)
  names(z) <- "volume"
  pred <- predFun(x)
  (y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - Q.025
}
c(uniroot(upper, interval = c(min(bladder$volume), x0.est), 
    tol = 1e-10, maxiter = 1000)$root,
  uniroot(lower, interval = c(x0.est, max(bladder$volume)), 
    tol = 1e-10, maxiter = 1000)$root)


## Results #####################################################################

# Wald-based interval:            (54.94, 131.89)
# Inversion interval:             (58.24, 137.05)
# PB percentile interval:         (58.37, 136.02)
# PB adjusted inversion interval: (57.94, 138.16)


## Figures #####################################################################

## Figure 1 (scatterplot)
setEPS()
postscript(paste(home.dir, "scatter.eps", sep="/"), width = 4, height = 4)
xyplot(HD ~ volume, groups = subject, data = bladder, 
       type = "b", lty = 2, #cex = 0.5, 
       xlab = "Volume (ml)",
       ylab = expression(paste("Height (mm) " %*% " depth (mm)")),
       scales = list(tck = c(1, 0)))
dev.off()

## Figure 2 (fitted mean response)
setEPS()
postscript(paste(home.dir, "fit.eps", sep="/"), width = 4, height = 4)
newx <- list(volume = seq(from = min(bladder$volume), 
                          to = max(bladder$volume), 
                          length = 500))
p1 <- xyplot(HD ~ volume, groups = subject, data = bladder)
xyplot(HD ~ volume, groups = subject, data = bladder, 
  type = "b", col = "lightgrey",
  xlab = "Volume (ml)", 
  ylab = expression(paste("Height (mm) " %*% " depth (mm)")),
  scales = list(tck = c(1, 0)), 
  panel = function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.lines(unlist(newx), predict(bladder.nlme, newx, level = 0), 
                col = "black", lwd = 2, lty = 2)
    panel.arrows(p1$x.limits[1], 70, x0.est, 70, angle = 25)   
    panel.arrows(x0.est, 70, x0.est, p1$y.limits[1], angle = 25)
})
# p3 <- plot(augPred(bladder.nlme, length.out = 101, level = 0:1), 
#   xlab = "Volume (ml)", 
#   ylab = expression(paste("Height (mm) " %*% " depth (mm)")))
# print(p2, pos = c(0, 0, 1/2, 1), more = TRUE) 
# print(p3, pos = c(1/2, 0, 1, 1), more = FALSE)
dev.off()

## Figure 3 (historgram of bootstrapped predictive pivot)
Q.boot <- as.numeric(bladder.pb$t[, 2])
setEPS()
postscript(paste(home.dir, "inversion-hist.eps", sep="/"), width = 7, 
           height = 3.5)
par(mfrow = c(1, 2), las = 1, cex.axis = 0.9)
hist(Q.boot, br = 50, freq = FALSE, col = "darkgrey", border = "white",
     xlab = "", main = "")
qqnorm(Q.boot, xlab = "Theoretical quantile", ylab = "Sample quantile",
       main = "")
qqline(v0.boot)
dev.off()

## Figure 4 (historgram of bootstrapped point estimate)
v0.boot <- as.numeric(bladder.pb$t[, 1])
setEPS()
postscript(paste(home.dir, "v0-hist.eps", sep="/"), width = 7, height = 3.5)
par(mfrow = c(1, 2), las = 1, cex.axis = 0.9)
hist(v0.boot, br = 50, freq = FALSE, col = "darkgrey", border = "white",
     xlab = "", main = "")
qqnorm(v0.boot, xlab = "Theoretical quantile", ylab = "Sample quantile",
       main = "")
qqline(v0.boot)
dev.off()

