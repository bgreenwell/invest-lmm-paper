################################################################################
# Setup
################################################################################

# Required packages
library(investr)
library(lme4)
library(nlme)
library(plyr)

# Load bladder volume data
load("/home/w108bmg/Desktop/Dropbox/Datasets/bladder.RData")

# Transform response
bladder$HD2 <- bladder$HD^(3/2)


################################################################################
# Fitted models
################################################################################

# Simple linear regression
lm_fit <- lm(HD2 ~ volume, data = bladder)

# Random intercept and slope model (nlme)
nlme_fit <- lme(HD2 ~ volume, random = list(subject = pdDiag( ~ volume)), 
               data = bladder)  # using nlme package

# Random intercept and slope model (lme4)
lme4_fit <- lmer(HD2 ~ volume + (0 + 1 | subject) + (0 + volume | subject), 
                data = bladder)  # using lme4 package

## Using invest function. Note that the CI for x0 based on the LMM is slightly
## shorter. This is not surprising since Var(Y0) = 17572.35 at x0.est for the 
## LMM, whereas, Var(Y0) is a constant 20249.88 in the LM.
x0.est <- invest(nlme_fit, y0 = 500, interval = "none")
invest(lm_fit, y0 = 500)
invest(nlme_fit, y0 = 500)
invest(lm_fit, y0 = 500, interval = "Wald")
invest(nlme_fit, y0 = 500, interval = "Wald")


################################################################################
# Simulate data
################################################################################

inrange <- function(x, lower, upper, strict = TRUE) {
  stopifnot(is.numeric(x) || is.numeric(lower) || is.numeric(upper) || 
              is.logical(strict))
  if (strict) {
    x > lower & x < upper
  } else {
    x >= lower & x <= upper
  }
}

# Calculate response variance at a particular predictor value
var_Y <- function(object, x) {
  vc <- as.data.frame(VarCorr(object))$sdcor
  vc[1]^2 + vc[2]^2*x^2 + vc[3]^2
}
sd_y0 <- sqrt(var_Y(lme4_fit, x = 8.015521))

# Function to generate data of different sample sizes
gen_data <- function(m, n) {
  volume <- rep(seq(from = 1, to = 17.5, length = n), times = m)
  subject <- rep(1:m, each = n)
  theta <- as.data.frame(VarCorr(lme4_fit))$sdcor
  beta <- getME(lme4_fit, "beta")
  alpha0 <- rnorm(m, mean = beta[1], sd = theta[1])
  alpha1 <- rnorm(m, mean = beta[2], sd = theta[2])
  HD <- rnorm(m * n, mean = alpha0[subject] + alpha1[subject] * volume,
              sd = theta[3])
  data.frame(subject, HD, volume)
}

# True unknown
true_x0 <- invest(nlme_fit, y0 = 500, interval = "none", mean.response = TRUE)

# Function to simulate confidence intervals for the true unknown
simulate_intervals <- function(m, n, nsim = 10000L) {
  raply(nsim, {
    d <- gen_data(30, 30)
    fit <- lme(HD ~ volume, random = list(subject = pdDiag( ~ volume)), data = d) 
    Y0 <- 500 + rnorm(1, mean = 0, sd = sd_y0)
    out <- tryCatch(invest(fit, y0 = Y0, lower = -10, upper = 50),
                    error = function(e) {
                      list("estimate" = NA, "lower" = NA, "upper" = NA)
                    })
    c("est" = out$estimate, "lwr" = out$lower, "upr" = out$upper)
  }, .progress = "text")
}


################################################################################
#
################################################################################

sim <- simulate_intervals(m = 30,n = 20, nsim = 10000L)
mean(inrange(true_x0, sim[, "lwr"], sim[, "upr"]))
