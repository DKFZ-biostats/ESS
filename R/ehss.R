ess <- function(prior,...) UseMethod("ess")
ess.default <- function(prior, ...) warning("ess does not know how to handle object")
#print.ehss<-function(x, ...) print(c(ecss=x$ecss), ...)

ess.normMix=function(prior, method = c("mix.moment","moment", "morita"), ...){
  if (is.null(sigma(prior))) stop("reference scale sigma has to be given in mixture object")
  method <- match.arg(method)
  if (method=="morita"|method=="moment") res=RBesT:::ess.normMix(prior,method,...)
  if (method=="mix.moment") res=sum(prior["w",]*sigma(prior)^2/prior["s",]^2)
  c(ess=res)
}


ess.betaMix=function(prior, method = c("mix.moment","moment", "morita"), ...){
  method <- match.arg(method)
  if (method=="morita"|method=="moment") res=RBesT:::ess.betaMix(prior,method,...)
  if (method=="mix.moment") res=sum(prior["w",]*(prior["a",]+prior["b",]))
  c(ess=res)
}



ehss <- function(prior,...) UseMethod("ehss")
ehss.default <- function(prior, ...) warning("ehss does not know how to handle object")
#print.ehss<-function(x, ...) print(c(ecss=x$ecss), ...)

ehss.normMix=function(prior, data,  n, m, se,  method = c("mix.moment","moment", "morita"), ...){
  if (!missing(data)) {
    m <- mean(data)
    n <- length(data)
    se <- sd(data)/sqrt(n)
  }
  else {
    if (missing(m) & (missing(n) | missing(se))) 
        stop("Either raw data or summary data (m and se) must be given.")
    if (!missing(se) & !missing(n)) {
        sigma(prior) <- se * sqrt(n)
        message(paste0("Updating reference scale to ", sigma(prior), 
            ".\nIt is recommended to use the sigma command instead.\nSee ?sigma or ?mixnorm."))
    }
    if (missing(se) & !missing(n)) {
        message("Using default prior reference scale ", sigma(prior))
        se <- sigma(prior)/sqrt(n)
    }
    if (!missing(se) & missing(n)) {
        message("Using default prior reference scale ", sigma(prior))
        n<- (sigma(prior)/se)^2
    }
  }
  
  # if (inherits(prior,"powerprior")) {
  #       pp=powerprior(prior, m=m, n=n, sigma=se * sqrt(n))
  #       res=ess(pp, method=method)
  # 
  # } 
  
  res=ess(postmix(prior, m=m, n=n, sigma=se * sqrt(n)), method=method, ...)-n

  c(ehss=res)
}


ehss.betaMix=function(prior, data, n, r, method = c("mix.moment","moment", "morita"), ...){
    if (!missing(data)) {
          assert_that(all(data %in% c(0, 1)))
          r <- sum(data)
          n <- length(data)
      }
  # if (inherits(prior,"powerprior")) {
  #       pp=powerprior(prior, r=r, n=n)
  #       res=ess(pp, method=method)
  # 
  # } else 
    
    res=ess(postmix(prior, r=r, n=n), method=method, ...)-n

    c(ehss=res)
}
