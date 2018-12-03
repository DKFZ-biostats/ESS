as.powerprior=function(prior, p.prior.a=1, p.prior.b=1){
  if (ncol(prior)!=1) stop("mixture object may only have one component")
  if (inherits(prior,"betaMix")) attr(prior,"p.prior")=c(p.prior.a,p.prior.b)
  class(prior)=c("powerprior",class(prior))
  prior
}



#############################
# powerparameter
  powerparameter<- function(prior,...) UseMethod("powerparameter")
  powerparameter.default <- function(prior, ...) warning("powerparameter does not know how to handle object")
  
  powerparameter.betaMix <- Vectorize(function(prior, n, r, p.prior.a=1, p.prior.b=1,...) {
    if (ncol(prior)!=1) stop("mixture object may only have one component")
    lik.d <- function(d) VGAM::dbetabinom.ab(r, n, p.prior.a + 
                                               sum(d * prior["a",]), p.prior.b + sum(d * prior["b",])) #VGAM::dbetabinom.ab
    opd <- my.spg(par = 0.005, #prior guess of d
                   fn = lik.d, 
                   lower = 0, upper = 1, 
                   control = list(maximize = TRUE, trace = FALSE,
                                  checkGrad=FALSE)) # BB::spg
    if (opd$convergence != 0) 
      print(opd)
    d <- opd$par
    return(d)
  },vectorize.args = "r")
  
  
  powerparameter.normMix <- Vectorize(function(prior, n, m, sigma,...) {  
    if (ncol(prior)!=1) stop("mixture object may only have one component")
    #m data mean
    # n data sample size
    # SD2 data variance
    #prior is a RBeST mix object
      
    var.mle=sigma^2/n #new mean variance
    var.prior=prior["s",]^2
    mean.prior=prior["m",]
    mle=m
    d <- var.prior/(max((mle-mean.prior)^2,var.mle+var.prior)-var.mle)
    return(d)
  },vectorize.args = "m")
  




####################
# powerprior
  powerprior <- function(prior,...) UseMethod("powerprior")
  powerprior.default <- function(prior, ...) warning("powerprior does not know how to handle object")
  
  powerprior.normMix=function(prior, n, m, sigma, ...){
    if (ncol(prior)!=1) stop("mixture object may only have one component")
  #  ds=ds_normal(X=X,N=N,SD2=sigma^2,prior)
    ds=powerparameter(prior,m=m,n=n,sigma=sigma)
    mixnorm(c(1,prior["m",],sqrt(prior["s",]^2/ds)),sigma=sigma)
  }
  # ehss can be calculated using ehss(powerprior(prior,m,n,sigma)) or ehss(postmix(as.powerprior(prior),m,n,sigma))
  
  powerprior.betaMix=function(prior, n, r, p.prior.a, p.prior.b, ...){
    if (missing(p.prior.a)| missing(p.prior.b)){
      if (!is.null(attr(prior,"p.prior"))){
        p.prior.a=attr(prior,"p.prior")[1]
        p.prior.b=attr(prior,"p.prior")[2]
      } else{
        p.prior.a=p.prior.b=1
      }
    }
    if (ncol(prior)!=1) stop("mixture object may only have one component")
  #  ds=ds_binomial(X=X,N=N,x=prior["a",1],n=prior["a",1]+prior["b",1])
    ds=powerparameter(prior,r=r,n=n,p.prior.a=p.prior.a,p.prior.b=p.prior.b)
    mixbeta(c(1,p.prior.a+ds*prior["a",],p.prior.b+ds*prior["b",]))
  }

  
  
################
# postmix
  postmix.powerprior=function(priormix, data, n, r, m, se, p.prior.a, p.prior.b, ...){
    if (inherits(priormix,"normMix")){
      if (!missing(data)) {
          m <- mean(data)
          n <- length(data)
          se <- sd(data)/sqrt(n)
      }
      else {
          if (missing(m) & (missing(n) | missing(se))) 
              stop("Either raw data or summary data (m and se) must be given.")
          if (!missing(se) & !missing(n)) {
              sigma(priormix) <- se * sqrt(n)
              message(paste0("Updating reference scale to ", sigma(priormix), 
                  ".\nIt is recommended to use the sigma command instead.\nSee ?sigma or ?mixnorm."))
          }
          if (missing(se) & !missing(n)) {
              message("Using default prior reference scale ", sigma(priormix))
              se <- sigma(priormix)/sqrt(n)
          }
          if (!missing(se) & missing(n)) {
              message("Using default prior reference scale ", sigma(priormix))
              n<- (sigma(priormix)/se)^2
          }
      }
      pp=powerprior(priormix, m=m, n=n, sigma=se * sqrt(n))
      return(postmix(pp,m=m,se=se))
    }
    if (inherits(priormix,"betaMix")){
      if (!missing(data)) {
          assert_that(all(data %in% c(0, 1)))
          r <- sum(data)
          n <- length(data)
      }
      if (missing(p.prior.a)| missing(p.prior.b)){
        if (!is.null(attr(priormix,"p.prior"))){
          p.prior.a=attr(priormix,"p.prior")[1]
          p.prior.b=attr(priormix,"p.prior")[2]
        } else{
          p.prior.a=p.prior.b=1
        }
      }
  
      pp=powerprior(priormix,r=r,n=n,p.prior.a=p.prior.a,p.prior.b=p.prior.b)
      return(postmix(pp,r=r,n=n))
    }
  }





# # legacy - will be deprecated .Deprecated
#   ds_binomial <-
#   Vectorize(function(X,N,x,n,p.prior.a=.001,p.prior.b=.001) {
#     lik.d <- function(d) VGAM::dbetabinom.ab(X, N, p.prior.a + 
#                                                sum(d * x), p.prior.b + sum(d * (n - x))) #VGAM::dbetabinom.ab
#     opd <- BB::spg(par = rep(0.005, length(x)), fn = lik.d, 
#                    lower = rep(0, length(x)), upper = rep(1, length(x)), 
#                    control = list(maximize = TRUE, trace = FALSE)) # BB::spg
#     if (opd$convergence != 0) 
#       print(opd)
#     d <- opd$par
#     return(d)
#   },vectorize.args = "X")
# 
# 
# # legacy - will be deprecated .Deprecated
#   ds_normal <- Vectorize(function(X,N,SD2,prior) {  # 
#     #X data mean
#     # N data sample size
#     # SD2 data variance
#     #prior is a RBeST mix object
#       
#     var.mle=SD2/N #new mean variance
#     var.prior=prior["s",]^2
#     mean.prior=prior["m",]
#     mle=X
#     d <- var.prior/(max((mle-mean.prior)^2,var.mle+var.prior)-var.mle)
#     return(d)
#   },vectorize.args = "X")



my.spg=function (par, fn, gr = NULL, method = 3, lower = -Inf, upper = Inf, 
  project = NULL, projectArgs = NULL, control = list(), quiet = FALSE, 
  alertConvergence = TRUE, ...) 
{
  box <- if (any(is.finite(upper))) 
    TRUE
  else if (any(is.finite(lower))) 
    TRUE
  else FALSE
  prj <- if (box) 
    TRUE
  else if (!is.null(project)) 
    TRUE
  else FALSE
  if (is.character(project)) 
    project <- get(project, mode = "function")
  if (box) {
    if (is.null(project)) {
      if (is.null(projectArgs)) 
        projectArgs <- list()
      if ((!is.null(projectArgs$lower)) | (!is.null(projectArgs$upper))) 
        warning("Using lower and upper spg arguments, ", 
          "not using those specified in projectArgs.")
      projectArgs$lower <- if (length(lower) == 1) 
        rep(lower, length(par))
      else lower
      projectArgs$upper <- if (length(upper) == 1) 
        rep(upper, length(par))
      else upper
      project <- function(par, lower, upper) {
        par[par < lower] <- lower[par < lower]
        par[par > upper] <- upper[par > upper]
        return(par)
      }
    }
    if (identical(project, BB::projectLinear)) {
      if ((!is.null(projectArgs$lower)) | (!is.null(projectArgs$upper))) 
        warning("Using lower and upper spg arguments, ", 
          "not using those specified in projectArgs.")
      if (is.null(projectArgs$A)) 
        stop("projectLinear requires the A matrix to be specified in projectArgs.")
      if (is.null(projectArgs$b)) 
        stop("projectLinear requires the b vector to be specified in projectArgs.")
      if (length(lower) == 1) 
        lower <- rep(lower, length(par))
      if (any(zi <- is.finite(lower))) {
        projectArgs$A <- rbind(projectArgs$A, diag(length(par))[zi, 
          ])
        projectArgs$b <- c(projectArgs$b, lower[zi])
      }
      if (length(upper) == 1) 
        upper <- rep(upper, length(par))
      if (any(zi <- is.finite(upper))) {
        projectArgs$A <- rbind(projectArgs$A, diag(-1, 
          length(par))[zi, ])
        projectArgs$b <- c(projectArgs$b, -upper[zi])
      }
    }
  }
  ctrl <- list(M = 10, maxit = 1500, ftol = 1e-10, gtol = 1e-05, 
    maxfeval = 10000, maximize = FALSE, trace = TRUE, triter = 10, 
    quiet = FALSE, eps = 1e-07, checkGrad = NULL, checkGrad.tol = 1e-06)
  namc <- names(control)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% 
      names(ctrl))])
  ctrl[namc] <- control
  M <- ctrl$M
  maxit <- ctrl$maxit
  gtol <- ctrl$gtol
  ftol <- ctrl$ftol
  maxfeval <- ctrl$maxfeval
  maximize <- ctrl$maximize
  trace <- ctrl$trace
  triter <- ctrl$triter
  eps <- ctrl$eps
  checkGrad <- ctrl$checkGrad
  checkGrad.tol <- ctrl$checkGrad.tol
  grNULL <- is.null(gr)
  fargs <- list(...)
  func <- if (maximize) 
    function(par, ...) c(-fn(par, ...))
  else function(par, ...) c(fn(par, ...))
  feval <- 1
  if (!grNULL){ # this line is the only differnce to BB::spg
    f.time <- system.time(f <- try(func(par, ...), silent = TRUE))
    if (is.null(checkGrad)) {
      if (((f.time[1] + f.time[2]) * 6 * length(par)) < 10) {
        checkGrad <- TRUE
      } else {
        checkGrad <- FALSE
        if (!grNULL) 
          warning("Default checking of gradient turned off because of time require.", 
            "See the help for spg to enable this.")
      }
    }
  } else f <- try(func(par, ...), silent = TRUE)
  if (class(f) == "try-error") 
    stop("Failure in initial function evaluation!", f)
  else if (!is.numeric(f) || 1 != length(f)) 
    stop("function must return a scalar numeric value!")
  else if (is.nan(f) | is.infinite(f) | is.na(f)) 
    stop("Failure in initial function evaluation!")
  f0 <- fbest <- f
  nmls <- function(p, f, d, gtd, lastfv, feval, func, maxfeval, 
    fargs) {
    gamma <- 1e-04
    fmax <- max(lastfv)
    alpha <- 1
    pnew <- p + alpha * d
    fnew <- try(do.call(func, append(list(pnew), fargs)), 
      silent = TRUE)
    feval <- feval + 1
    if (class(fnew) == "try-error" | is.nan(fnew)) 
      return(list(p = NA, f = NA, feval = NA, lsflag = 1))
    while (fnew > fmax + gamma * alpha * gtd) {
      if (alpha <= 0.1) 
        alpha <- alpha/2
      else {
        atemp <- -(gtd * alpha^2)/(2 * (fnew - f - alpha * 
          gtd))
        if (atemp < 0.1 | atemp > 0.9 * alpha) 
          atemp <- alpha/2
        alpha <- atemp
      }
      pnew <- p + alpha * d
      fnew <- try(do.call(func, append(list(pnew), fargs)), 
        silent = TRUE)
      feval <- feval + 1
      if (class(fnew) == "try-error" | is.nan(fnew)) 
        return(list(p = NA, f = NA, feval = NA, lsflag = 1))
      if (feval > maxfeval) 
        return(list(p = NA, f = NA, feval = NA, lsflag = 2))
    }
    return(list(p = pnew, f = fnew, feval = feval, lsflag = 0))
  }
  if (!grNULL & checkGrad) {
    requireNamespace("numDeriv", quietly = TRUE)
    grad.num <-numDeriv::grad(x = par, func = fn, ...)
    grad.analytic <- gr(par, ...)
    max.diff <- max(abs((grad.analytic - grad.num)/(1 + 
      abs(fn(par, ...)))))
    if (!max.diff < checkGrad.tol) {
      cat("Gradient check details:  max. relative difference in gradients= ", 
        max.diff, "\n\n  analytic gradient:", grad.analytic, 
        "\n\n  numerical gradient:", grad.num)
      stop("Analytic gradient does not seem correct! See comparison above. ", 
        "Fix it, remove it, or increase checkGrad.tol.")
    }
  }
  if (grNULL) 
    gr <- function(par, ...) {
      df <- rep(NA, length(par))
      for (i in 1:length(par)) {
        dx <- par
        dx[i] <- dx[i] + eps
        df[i] <- (func(dx, ...) - f)/eps
      }
      df
    }
  lmin <- 1e-30
  lmax <- 1e+30
  iter <- 0
  lastfv <- rep(-1e+99, M)
  fbest <- NA
  fchg <- Inf
  grad <- if (maximize & !grNULL) 
    function(par, ...) -gr(par, ...)
  else function(par, ...) gr(par, ...)
  if (prj) {
    par <- try(do.call(project, append(list(par), projectArgs)), 
      silent = TRUE)
    if (class(par) == "try-error") 
      stop("Failure in projecting initial guess!", par)
  }
  if (any(is.nan(par), is.na(par))) 
    stop("Failure in initial guess!")
  pbest <- par
  g <- try(grad(par, ...), silent = TRUE)
  if (class(g) == "try-error") 
    stop("Failure in initial gradient evaluation!", g)
  else if (any(is.nan(g))) 
    stop("Failure in initial gradient evaluation!")
  lastfv[1] <- fbest <- f
  pg <- par - g
  if (prj) {
    pg <- try(do.call(project, append(list(pg), projectArgs)), 
      silent = TRUE)
    if (class(pg) == "try-error") 
      stop("Failure in initial projection!", pg)
  }
  if (any(is.nan(pg))) 
    stop("Failure in initial projection!")
  pg <- pg - par
  pg2n <- sqrt(sum(pg * pg))
  pginfn <- max(abs(pg))
  gbest <- pg2n
  if (pginfn != 0) 
    lambda <- min(lmax, max(lmin, 1/pginfn))
  if (trace) 
    cat("iter: ", 0, " f-value: ", f0 * (-1)^maximize, " pgrad: ", 
      pginfn, "\n")
  lsflag <- NULL
  while (pginfn > gtol & iter <= maxit & fchg > ftol) {
    iter <- iter + 1
    d <- par - lambda * g
    if (prj) {
      d <- try(do.call(project, append(list(d), projectArgs)), 
        silent = TRUE)
      if (class(d) == "try-error" | any(is.nan(d))) {
        lsflag <- 4
        break
      }
    }
    d <- d - par
    gtd <- sum(g * d)
    if (is.infinite(gtd)) {
      lsflag <- 4
      break
    }
    nmls.ans <- nmls(par, f, d, gtd, lastfv, feval, func, 
      maxfeval, fargs)
    lsflag <- nmls.ans$lsflag
    if (lsflag != 0) 
      break
    fchg <- abs(f - nmls.ans$f)
    f <- nmls.ans$f
    pnew <- nmls.ans$p
    feval <- nmls.ans$feval
    lastfv[(iter%%M) + 1] <- f
    gnew <- try(grad(pnew, ...), silent = TRUE)
    if (class(gnew) == "try-error" | any(is.nan(gnew))) {
      lsflag <- 3
      break
    }
    s <- pnew - par
    y <- gnew - g
    sts <- sum(s * s)
    yty <- sum(y * y)
    sty <- sum(s * y)
    if (method == 1) 
      lambda <- if (sts == 0 | sty < 0) 
        lmax
      else min(lmax, max(lmin, sts/sty))
    else if (method == 2) 
      lambda <- if (sty < 0 | yty == 0) 
        lmax
      else min(lmax, max(lmin, sty/yty))
    else if (method == 3) 
      lambda <- if (sts == 0 | yty == 0) 
        lmax
      else min(lmax, max(lmin, sqrt(sts/yty)))
    par <- pnew
    g <- gnew
    pg <- par - g
    if (prj) {
      pg <- try(do.call(project, append(list(pg), projectArgs)), 
        silent = TRUE)
      if (class(pg) == "try-error" | any(is.nan(pg))) {
        lsflag <- 4
        break
      }
    }
    pg <- pg - par
    pg2n <- sqrt(sum(pg * pg))
    pginfn <- max(abs(pg))
    f.rep <- (-1)^maximize * f
    if (trace && (iter%%triter == 0)) 
      cat("iter: ", iter, " f-value: ", f.rep, " pgrad: ", 
        pginfn, "\n")
    if (f < fbest) {
      fbest <- f
      pbest <- pnew
      gbest <- pginfn
    }
  }
  if (is.null(lsflag)) {
    if (!quiet) 
      warning("convergence tolerance satisified at intial parameter values.")
    lsflag <- 0
  }
  if (lsflag == 0) {
    if (pginfn <= gtol | fchg <= ftol) 
      conv <- list(type = 0, message = "Successful convergence")
    if (iter >= maxit) 
      conv <- list(type = 1, message = "Maximum number of iterations exceeded")
    f.rep <- (-1)^maximize * fbest
    par <- pbest
  }
  else {
    par <- pbest
    f.rep <- f <- (-1)^maximize * fbest
    pginfn <- gbest
    if (lsflag == 1) 
      conv <- list(type = 3, message = "Failure:  Error in function evaluation")
    if (lsflag == 2) 
      conv <- list(type = 2, message = "Maximum function evals exceeded")
    if (lsflag == 3) 
      conv <- list(type = 4, message = "Failure:  Error in gradient evaluation")
    if (lsflag == 4) 
      conv <- list(type = 5, message = "Failure:  Error in projection")
  }
  if (alertConvergence && (0 != conv$type)) 
    warning("Unsuccessful convergence.")
  return(list(par = par, value = f.rep, gradient = pginfn, 
    fn.reduction = (-1)^maximize * (f0 - f), iter = iter, 
    feval = feval, convergence = conv$type, message = conv$message))
}



