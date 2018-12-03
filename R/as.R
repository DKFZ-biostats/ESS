as.mix <- function(prior,...) UseMethod("as.mix")
as.mix.default <- function(prior, ...) warning("as.mix does not know how to handle object")

as.mix.beta=function(prior,...){
  do.call(mixbeta,as.list(data.frame(t(cbind(attr(prior,"weights"),attr(prior,"pars"))))))
}
as.mix.normal=function(prior,...){
  do.call(mixnorm,as.list(data.frame(t(cbind(attr(prior,"weights"),attr(prior,"pars"))))))
}



as.mixture.prior <- function(prior,...) UseMethod("as.mixture.prior")
as.mixture.prior.default <- function(prior, ...) warning("as.mixture.prior does not know how to handle object")

as.mixture.prior.normMix=function(prior,...){
  requireNamespace("StudyPrior", quietly = TRUE)
  StudyPrior::create.mixture.prior(type="normal", pars=t(prior[-1,]), weights = prior[1,])
}
as.mixture.prior.betaMix=function(prior,...){
 requireNamespace("StudyPrior", quietly = TRUE)
  StudyPrior::create.mixture.prior(type="beta", pars=t(prior[-1,]), weights = prior[1,])
}




