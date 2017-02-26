## Functions -------------------------------------------------------------------

## bootMer2
bootMer2 <- function (x, FUN, FUN0, nsim = 1, seed = NULL, use.u = FALSE, 
                      type = c("parametric", "semiparametric"), verbose = FALSE, 
                      .progress = "none", PBargs = list(), 
                      parallel = c("no", "multicore", "snow"), 
                      ncpus = getOption("boot.ncpus", 1L), cl = NULL) 
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (.progress != "none") {
    pbfun <- get(paste0(.progress, "ProgressBar"))
    setpbfun <- get(paste0("set", .simpleCap(.progress), 
                           "ProgressBar"))
    pb <- do.call(pbfun, PBargs)
  }
  if (missing(parallel)) 
    parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
  }
  do_parallel <- (ncpus > 1L && (have_mc || have_snow))
  if (do_parallel & .progress != "none") 
    message("progress bar disabled for parallel operations")
  FUN <- match.fun(FUN)
  type <- match.arg(type)
  if (!is.null(seed)) 
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  mc <- match.call()
  t0 <- FUN0(x)
  if (!is.numeric(t0)) 
    stop("bootMer currently only handles functions that return numeric vectors")
  mle <- list(beta = getME(x, "beta"), theta = getME(x, "theta"))
  if (isLMM(x)) 
    mle <- c(mle, list(sigma = sigma(x)))
  if (type == "parametric") {
    ss <- simulate(x, nsim = nsim, use.u = use.u, na.action = na.exclude)
  }
  else {
    if (use.u) {
      if (isGLMM(x)) 
        warning("semiparametric bootstrapping is questionable for GLMMs")
      ss <- replicate(nsim, fitted(x) + sample(residuals(x, 
                                                         "response")), simplify = FALSE)
    }
    else {
      stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
    }
  }
  ffun <- local({
    FUN
    refit
    x
    ss
    verbose
    do_parallel
    length.t0 <- length(t0)
    function(i) {
      foo <- try(FUN(refit(x, ss[[i]])), silent = TRUE)
      if (verbose) {
        cat(sprintf("%5d :", i))
        str(foo)
      }
      if (!do_parallel && .progress != "none") {
        setpbfun(pb, i/nsim)
      }
      if (inherits(foo, "try-error")) 
        rep(NA, length.t0)
      else foo
    }
  })
  simvec <- seq_len(nsim)
  res <- if (do_parallel) {
    if (have_mc) {
      parallel::mclapply(simvec, ffun, mc.cores = ncpus)
    }
    else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", 
                                             ncpus))
        parallel::clusterExport(cl, varlist = getNamespaceExports("lme4"))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, simvec, ffun)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, simvec, ffun)
    }
  }
  else lapply(simvec, ffun)
  t.star <- do.call(cbind, res)
  rownames(t.star) <- names(t0)
  if ((numFail <- sum(apply(is.na(t.star), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", 
            nsim, ")")
  }
  s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame, 
                      seed = .Random.seed, statistic = FUN, sim = "parametric", 
                      call = mc, ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle), 
                 class = "boot")
  attr(s, "bootFail") <- numFail
  s
}

## This is a safer version of the uniroot function suitable for simulations in
## which the position of the root is random. If at first the root is not 
## contained in the interval (a, b), then the interval is extended in one or 
## both directions until this condition is met or maxk expansions have been
## performed.
uniroot2 <- function (f, interval, ..., lower = min(interval), 
                            upper = max(interval), f.lower = f(lower, ...), 
                            f.upper = f(upper, ...), 
                            tol = .Machine$double.eps^0.25, maxiter = 1000,
                            extend = c("none", "left", "right", "both"), 
                            frac = 0.05, maxk = 1000) 
{  
  extend <- match.arg(extend)
  if (!missing(interval) && length(interval) != 2L) {
    stop("'interval' must be a vector of length 2")
  }
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper) {
    stop("lower < upper  is not fulfilled")
  }
  if (is.na(f.lower)) {
    stop("f.lower = f(lower) is NA")
  }
  if (is.na(f.upper)) {
    stop("f.upper = f(upper) is NA")
  }
  if (f.lower * f.upper > 0) {
    if (extend != "none") {
      k <- 1
      if (extend == "left") { # extend lower bound
        repeat {
          lower <- lower - frac
          f.lower = f(lower, ...)
          k <- k + 1
          if (f.lower * f.upper <= 0 || k >= maxk) break
        }
      } else if (extend == "right") { # extend upper bound
        repeat {
          upper <- upper + frac
          f.upper = f(upper, ...)
          k <- k + 1
          if (f.lower * f.upper <= 0 || k >= maxk) break
        }
      } else { # extend lower bound and upper bound
        repeat {
          lower <- lower - frac
          upper <- upper + frac
          f.lower = f(lower, ...)
          f.upper = f(upper, ...)
          k <- k + 1
          if (f.lower * f.upper <= 0 || k >= maxk) break
        }
      }
    } else {
      stop("f() values at end points not of opposite sign")
    }
  }
  uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = maxiter)
}


## Cube root function
curt <- function(x) sign(x) * abs(x)^(1/3)

## Function to generate data
genData <- function(n = params$n, m = params$m, beta = params$beta, 
                    theta = params$theta) {
  subject <- rep(1:m, each = n)
  x <- rep(seq(from = 0, to = 1, length = n), times = m) 
  B0 <- rnorm(m, mean = 0, sd = sqrt(theta[1]))
  B1 <- rnorm(m, mean = 0, sd = sqrt(theta[2]))
  y <- rnorm(m*n, mean = beta[1]+B0[subject] + 
               (beta[2]+B1[subject])*x + beta[3]*curt(x), 
             sd = sqrt(theta[3]))
  data.frame(x = x, y = y, subject = factor(subject))
}

## Function to calculate inverse estimate given a mixed model object
xest <- function(object, y0 = params$y0, lower, upper, tol = 1e-10, 
                 maxiter = 1000, extend = "both", ...) {
  fixed <- as.numeric(fixef(object))
  f <- function(x) fixed[1] + fixed[2]*x + fixed[3]*curt(x) - y0
  uniroot2(f, lower = lower, upper = upper, tol = tol, 
                 maxiter = maxiter, extend = extend, ...)$root
}

## Function to summarize confidence intervals
mc.summary <- function(x, boot = FALSE) {
  
  ## Function that returns coverage and length of a single C.I.
  coverage.and.length <- function(x) {
    .cov <- if (x[1] <= params$x0 && params$x0 <= x[2]) 1 else 0
    .len <- x[2] - x[1]
    c(.cov, .len)
  }
  
  ## Summarize non bootstrap intervals
  if (!boot) {
    res <- apply(x, 1, coverage.and.length)
    res <- c(mean(res[1, ]), mean(res[2, ]), sd(res[2, ]))
    names(res) <- c("Coverage", "Length", "StDev(Length)")
    res
    
    ## Summarize bootstrap intervals
  } else {
    d <- ldply(x, function(x) {
      c(coverage.and.length(x[1, ]), # normal
        coverage.and.length(x[2, ]), # basic
        coverage.and.length(x[3, ]), # percentile
        coverage.and.length(x[4, ]), # adjusted inversion 11
        coverage.and.length(x[5, ]), # adjusted inversion 22
        coverage.and.length(x[6, ]), # adjusted inversion 12
        coverage.and.length(x[7, ])) # adjusted inversion 21   
    })
    boot.norm  <- c(mean(d[, 1]),  mean(d[, 2]),  sd(d[, 2]))
    boot.basic <- c(mean(d[, 3]),  mean(d[, 4]),  sd(d[, 4]))
    boot.perc  <- c(mean(d[, 5]),  mean(d[, 6]),  sd(d[, 6]))
    boot.inv11 <- c(mean(d[, 7]),  mean(d[, 8]),  sd(d[, 8]))
    boot.inv22 <- c(mean(d[, 9]),  mean(d[, 10]), sd(d[, 10]))
    boot.inv12 <- c(mean(d[, 11]), mean(d[, 12]), sd(d[, 12]))
    boot.inv21 <- c(mean(d[, 13]), mean(d[, 14]), sd(d[, 14]))
    res <- rbind(boot.norm, boot.basic, boot.perc, boot.inv11, boot.inv22,
                 boot.inv12, boot.inv21)
    colnames(res) <- c("Coverage", "Length", "StDev(Length)")
    res
  }
  
}

## Function to calculate the Wald-based C.I.
wald <- function(d, Y0, lower = params$x0-0.5, upper = params$x0+0.5) {
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), data = d) 
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  x0.est <- xest(mod, y0 = Y0, lower = params$x0-2, upper = params$x0+2)
  var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
    summary(mod)$sigma^2
  
  ## Inverse estimator as a function of the parameters
  dmFun <- function(pars) {
    fun <- function(x) (t(c(1, x, curt(x))) %*% pars[-length(pars)]) - 
      pars[length(pars)]
    uniroot2(fun, lower = lower, upper = upper, tol = 1e-10, 
                   maxiter = 1000, extend = "both")$root
  }
  
  ## Assign parameter names, calculate gradient, and return standard error
  covmat <- diag(4)
  covmat[1:3, 1:3] <- as.matrix(vcov(mod))
  covmat[4, 4] <- var.y0
  pars <- c(fixef(mod), Y0)
  gv <- attr(numericDeriv(quote(dmFun(pars)), "pars"), "gradient")
  se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
  x0.est + qnorm(c(0.025, 0.975))*se
  
}

## Function to calculate the inversion interval
inversion <- function(d, Y0, q1 = qnorm(0.025), q2 = qnorm(0.975), 
                      lower = params$x0-0.5, upper = params$x0+0.5, 
                      ...) {
  
  ## Fit model using nlme package and estimate x0 and Var(Y0)
  mod <- lme(y ~ x + curt(x), random = list(subject = pdDiag(~x)), data = d)   
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  x0.est <- xest(mod, y0 = Y0, lower = lower, upper = upper)
  var.y0 <- getVarCov(mod)[1, 1] + getVarCov(mod)[2, 2]*x0.est^2 + 
    summary(mod)$sigma^2
  
  ## FIXME: can the following code be fixed to avoid errors in the bootstrap
  ##        simulation? The problem seems to occur in finding the roots of fun1 
  ##        and fun2.
  
  ## Prediction function that also returns standard error
  predFun <- function(x) {
    z <- list("x" = x)
    fit <- predict(mod, newdata = z, level = 0)
    se.fit <- sqrt(diag(cbind(1, unlist(z), curt(unlist(z))) %*% mod$varFix %*% 
                          t(cbind(1, unlist(z), curt(unlist(z))))))
    list(fit = fit, se.fit = se.fit)
  }
  
  ## Inverse functions for calculating confidence limits
  fun1 <- function(x) { 
    pred <- predFun(x)
    (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q2
  }
  fun2 <- function(x) { 
    pred <- predFun(x)
    (Y0 - pred$fit)/sqrt((var.y0 + pred$se.fit^2)) - q1
  }
  
  ## Check curves
  curve(fun1, from = x0.est-5, to = x0.est+10, lwd = 2, ylab = "Bounds", 
        ylim = c(-10, 10))
  curve(fun2, lwd = 2, col = "red", add = TRUE)
  abline(h = 0, v = x0.est, lty = 2)
  
  ## Find roots based on newer bounds
  .lower <- try(uniroot2(fun1, interval = c(lower, x0.est), tol = 1e-10, 
                               maxiter = 1000, extend = "left", ...)$root,
                silent = TRUE)
  .upper <- try(uniroot2(fun2, interval = c(x0.est, upper), tol = 1e-10, 
                               maxiter = 1000, extend = "right", ...)$root, 
                silent = TRUE)
  if (inherits(.lower, "try-error")) .lower <- -Inf
  if (inherits(.upper, "try-error")) .upper <- Inf
  
  ## Find roots of inverse function
  c(.lower, .upper)
  
}

## Function to calculate bootstrap intervals
parboot <- function(d, Y0, R = 999, .parallel = TRUE) {
  
  ## Fit model using lme4 package and estimate x0 and Var(Y0)
  mod <- lmer(y ~ x + curt(x) + (0+1|subject) + (0+x|subject), data = d)
  if (missing(Y0)) Y0 <- rnorm(1, mean = params$y0, sd = sqrt(params$var.y0))
  x0.est <- xest(mod, y0 = Y0, lower = params$x0-2, upper =  params$x0+2)
  var.y0 <- VarCorr(mod)[[1]][1] + VarCorr(mod)[[2]][1]*x0.est^2 + 
    sigma(mod)^2
  
  ## Function to calculate bootstrap estimates
  bootFun <- function(.) {
    
    ## Calculate bootstrap inverse estimate
    y0.boot <- rnorm(1, mean = Y0, sd = sqrt(var.y0))
    x0.boot <- xest(., y0 = y0.boot, lower = params$x0-2, upper =  params$x0+2)
    
    ## Calculate bootstrap predictive pivots
    covb <- as.matrix(vcov(.))
    beta.boot <- as.numeric(fixef(.))
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, curt(x0.est)))) 
    ## FIXME: Should variances be calculated at x0.est or x0.boot?
    var1.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var1.mu0 <- t(c(1, x0.est, curt(x0.est))) %*% covb %*% 
      c(1, x0.est, curt(x0.est))
    var2.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.boot^2 + sigma(.)^2
    var2.mu0 <- t(c(1, x0.boot, curt(x0.boot))) %*% covb %*% 
      c(1, x0.boot, curt(x0.boot))
    Q11.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var1.mu0)
    Q22.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var2.mu0)
    Q12.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var2.mu0)
    Q21.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var1.mu0)
    
    ## Return bootstrap estimates
    c(x0.boot, Q11.boot, Q22.boot, Q12.boot, Q21.boot)
    
  }
  
  ## Function to calculate bootstrap estimates
  bootFun0 <- function(.) {
    
    ## Calculate bootstrap inverse estimate
    y0.boot <- Y0 
    x0.boot <- xest(., y0 = y0.boot, lower = params$x0-2, upper =  params$x0+2)
    
    ## Calculate bootstrap predictive pivot
    covb <- as.matrix(vcov(.))        # (X' V^-1 X)^-1
    beta.boot <- as.numeric(fixef(.))
    mu0.boot <-  as.numeric(crossprod(beta.boot, c(1, x0.est, curt(x0.est))))
    ## FIXME: Should variances be calculated at x0.est or x0.boot?
    var1.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.est^2 + sigma(.)^2
    var1.mu0 <- t(c(1, x0.est, curt(x0.est))) %*% covb %*% 
      c(1, x0.est, curt(x0.est))
    var2.y0 <- VarCorr(.)[[1]][1] + VarCorr(.)[[2]][1]*x0.boot^2 + sigma(.)^2
    var2.mu0 <- t(c(1, x0.boot, curt(x0.boot))) %*% covb %*% 
      c(1, x0.boot, curt(x0.boot))
    Q11.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var1.mu0)
    Q22.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var2.mu0)
    Q12.boot <- (y0.boot - mu0.boot)/sqrt(var1.y0+var2.mu0)
    Q21.boot <- (y0.boot - mu0.boot)/sqrt(var2.y0+var1.mu0)
    
    ## Return bootstrap estimates
    c(x0.boot, Q11.boot, Q22.boot, Q12.boot, Q21.boot)
    
  }
  
  ## Generate bootstrap samples
  x0.pb <- if (.parallel) {
    bootMer2(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R, 
             parallel = "multicore", ncpus = 4)
  } else {
    bootMer2(mod, FUN = bootFun, FUN0 = bootFun0, nsim = R)
  }
  
  ## Calculate bootstrap CIs
  x0.pb.ci <- boot.ci(x0.pb, type = c("norm", "basic", "perc"))
  q11 <- unname(quantile(x0.pb$t[, 2], c(0.025, 0.975), na.rm = TRUE))
  q22 <- unname(quantile(x0.pb$t[, 3], c(0.025, 0.975), na.rm = TRUE))
  q12 <- unname(quantile(x0.pb$t[, 4], c(0.025, 0.975), na.rm = TRUE))
  q21 <- unname(quantile(x0.pb$t[, 5], c(0.025, 0.975), na.rm = TRUE))
  
  rbind(x0.pb.ci$normal[2:3],
        x0.pb.ci$basic[4:5],
        x0.pb.ci$perc[4:5],
        inversion(d, q1 = q11[1], q2 = q11[2], Y0 = Y0),
        inversion(d, q1 = q22[1], q2 = q22[2], Y0 = Y0),
        inversion(d, q1 = q12[1], q2 = q12[2], Y0 = Y0),
        inversion(d, q1 = q21[1], q2 = q21[2], Y0 = Y0))  
  
}
