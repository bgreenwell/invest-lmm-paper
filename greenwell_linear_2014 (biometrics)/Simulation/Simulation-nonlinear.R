## Clear working directory -----------------------------------------------------
rm(list = ls())

## Load packages and source code -----------------------------------------------
library(lattice) # for plotting grouped data
library(plyr)    # for working with the data
library(lme4)    # for fitting mixed models and bootstrap inference
library(nlme)    # for fitting mixed models
library(boot)    # for calculating bootstrap confidence intervals
source("/home/w108bmg/Desktop/Simulation/utility.R")

## Simulation parameters -------------------------------------------------------

## Parameters
params <- list(
  nsim = 1000,                  # simulation size
  m = 15,                       # number of subjects
  n = 10,                       # sample size per subject
  beta = c(0, 0.1, 0.9),        # fixed effecs
  theta = c(0.01, 0.05, 0.001), # variance components
  x0 = 1                     # true unknown
)
params$y0 <- params$beta[1] + params$beta[2]*params$x0 + 
  params$beta[3]*curt(params$x0)
params$var.y0 <- params$theta[1] + params$theta[2]*params$x0^2 + 
  params$theta[3]

## Sample data frame
set.seed(101)
simdata <- genData(n = 10, m = 15)
p <- xyplot(y ~ x, , data = simdata, type = "n")
xyplot(y ~ x, groups = subject, data = simdata, type = "b", 
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...) 
         panel.curve(params$beta[1] + params$beta[2]*x + params$beta[3]*curt(x), 
                     lwd = 3, from = 0, to = 1)
         panel.segments(p$x.limits[1], params$y0, params$x0, params$y0)
         panel.arrows(params$x0, params$y0, params$x0, p$y.limits[1]) 
       })

## Sample models
fit.nlme <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), 
               data = simdata)
fit.lme4 <- lmer(y ~ x + curt(x) + (0+1|subject) + (0+x|subject), 
                 data = simdata)

## Compare point estimates (both should be the same)
est1 <- xest(fit.nlme, lower = 0, upper = 1)
est2 <- xest(fit.lme4, lower = 0, upper = 1)
c(est1, est2)
all.equal(est1, est2)

## Compare confidence limits
wald(simdata, Y0 = params$y0)
inversion(simdata, Y0 = params$y0)
parboot(simdata, Y0 = params$y0, R = 9999)
# system.time(z1 <- parboot(simdata, Y0 = params$y0))
# system.time(z2 <- parboot(simdata, Y0 = params$y0, .parallel = FALSE))

## Simulations -----------------------------------------------------------------

## Path for savinf results
path <- "/home/w108bmg/Desktop/Simulation/Nonlinear results/1.00"

## Simulate data frames
set.seed(5746)
dfs <- rlply(params$nsim, genData)

## Simulation for the Wald-based interval --------------------------------------
mc.wald <- laply(dfs, wald, .progress = "text")
mc.summary(mc.wald)
save(mc.wald, file = paste(path, "mc.wald.RData", sep = "/"))

## Simulation for the inversion interval ---------------------------------------
mc.inversion <- laply(dfs, inversion, .progress = "text")
mc.summary(mc.inversion)
save(mc.inversion, file = paste(path, "mc.inversion.RData", sep = "/"))

## Simulation for the PB intervals (~ 8 hrs) -----------------------------------
# mc.parboot <- llply(dfs, parboot, .progress = "text")
# round(mc.summary(mc.parboot, boot = TRUE), 4)
# save(mc.parboot, file = paste(path, "mc.parboot.RData", sep = "/"))

## Splitting it up is safer!
mc.parboot1 <- llply(dfs[1:100], parboot, .progress = "text")
mc.summary(mc.parboot1, boot = TRUE)
save(mc.parboot1, file = paste(path, "mc.parboot1.RData", sep = "/"))
  
mc.parboot2 <- llply(dfs[101:200], parboot, .progress = "text")
mc.summary(mc.parboot2, boot = TRUE)
save(mc.parboot2, file = paste(path, "mc.parboot2.RData", sep = "/"))
     
mc.parboot3 <- llply(dfs[201:300], parboot, .progress = "text") 
mc.summary(mc.parboot3, boot = TRUE)
save(mc.parboot3, file = paste(path, "mc.parboot3.RData", sep = "/"))
     
mc.parboot4 <- llply(dfs[301:400], parboot, .progress = "text") 
mc.summary(mc.parboot4, boot = TRUE)
save(mc.parboot4, file = paste(path, "mc.parboot4.RData", sep = "/"))
     
mc.parboot5 <- llply(dfs[401:500], parboot, .progress = "text")
mc.summary(mc.parboot5, boot = TRUE)
save(mc.parboot5, file = paste(path, "mc.parboot5.RData", sep = "/"))
     
mc.parboot6 <- llply(dfs[501:600], parboot, .progress = "text") 
mc.summary(mc.parboot6, boot = TRUE)
save(mc.parboot6, file = paste(path, "mc.parboot6.RData", sep = "/"))
     
mc.parboot7 <- llply(dfs[601:700], parboot, .progress = "text")
mc.summary(mc.parboot7, boot = TRUE)
save(mc.parboot7, file = paste(path, "mc.parboot7.RData", sep = "/"))
     
mc.parboot8 <- llply(dfs[701:800], parboot, .progress = "text")
mc.summary(mc.parboot8, boot = TRUE)
save(mc.parboot8, file = paste(path, "mc.parboot8.RData", sep = "/"))
     
mc.parboot9 <- llply(dfs[801:900], parboot, .progress = "text")
mc.summary(mc.parboot9, boot = TRUE)
save(mc.parboot9, file = paste(path, "mc.parboot9.RData", sep = "/"))
     
mc.parboot10 <- llply(dfs[901:1000], parboot, .progress = "text")
mc.summary(mc.parboot10, boot = TRUE)
save(mc.parboot10, file = paste(path, "mc.parboot10.RData", sep = "/"))
     
mc.parboot <- c(mc.parboot1, mc.parboot2, mc.parboot3, mc.parboot4, mc.parboot5, 
                mc.parboot6, mc.parboot7, mc.parboot8, mc.parboot9, mc.parboot10)
mc.summary(mc.parboot, boot = TRUE)
save(mc.parboot, file = paste(path, "mc.parboot.RData", sep = "/"))
