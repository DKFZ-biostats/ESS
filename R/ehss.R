ehss <- function(prior,...) UseMethod("ehss")
ehss.default <- function(prior, ...) warning("ehss does not know how to handle object")
#print.ehss<-function(x, ...) print(c(ecss=x$ecss), ...)

ehss.normMix=function(prior, method = c("mix.moment","moment", "morita"), ...){
  if (is.null(sigma(prior))) stop("reference scale sigma has to be given in mixture object")
  method <- match.arg(method)
  if (method=="morita"|method=="moment") res=RBesT::ess(prior,method,...)
  if (method=="mix.moment") res=sum(prior["w",]*sigma(prior)^2/prior["s",]^2)
  c(ehss=res)
}


ehss.betaMix=function(prior, method = c("mix.moment","moment", "morita"), ...){
  method <- match.arg(method)
  if (method=="morita"|method=="moment") res=RBesT::ess(prior,method,...)
  if (method=="mix.moment") res=sum(prior["w",]*(prior["a",]+prior["b",]))
  c(ehss=res)
}
