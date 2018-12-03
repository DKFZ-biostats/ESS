freqError <- function(prior,...) UseMethod("freqError")
freqError.default <- function(prior, ...) warning("freqError does not know how to handle object")


freqError.betaMix=freqErrorMix.betaMix=function(prior,y,k,true.mean, treat, test, ...){
#  dec <- decision1S(0.975, treat, lower.tail=FALSE)
  posterior=postmix(prior,r=y,n=k) 
  abs((pmix(posterior,treat,lower.tail = T)>0.975)-test)
}


freqError.normMix=function(prior,y,k,true.mean, treat, test, ...){
  sigma2=RBesT::sigma(prior)^2
  var=  1/(1/prior["s",1]^2+k/sigma2)
  mean= var*(prior["m",1]/(prior["s",1]^2)+y/(sigma2/k))
  abs((pnorm(treat,mean=mean,sd=sqrt(var),lower.tail = T)>0.975)-test)

}



freqErrorMix.normMix=function(prior,y,k,true.mean, treat, test, ...){
#  posterior= my.postmix.normMix(prior, m=y,n=k,se=RBesT::sigma(prior)/sqrt(k))
#  posterior= lapply(y, function(yy) suppressMessages(suppressWarnings(postmix(prior, m=yy,n=k,se=RBesT::sigma(prior)/sqrt(k)))))
  posterior= lapply(y, function(yy) suppressMessages(suppressWarnings(my.postmix.normMix(prior, m=yy,n=k,se=RBesT::sigma(prior)/sqrt(k)))))
  posterior=lapply(posterior,function(post)
    my.mixnorm(post))
  
  res=sapply(posterior,function(post){
    post1=post
    abs((pmix(post1,treat,lower.tail = T)>0.975)-test)
  })
res

}

my.mixnorm=function (mat) {
    rownames(mix) <- c("w", "m", "s")
#    assert_that(all(mix["s", ] > 0))
    # if (!missing(sigma)) {
    #     assert_number(sigma, lower = 0)
    #     attr(mix, "sigma") <- sigma
    # }
    class(mix) <- c("normMix", "mix")
    likelihood(mix) <- "normal"
    mix
}
