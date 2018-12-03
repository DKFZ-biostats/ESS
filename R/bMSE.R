BMSE <- function(prior,...) UseMethod("BMSE")
BMSE.default <- function(prior, ...) warning("BMSE does not know how to handle object")


BMSE.mix <-function(prior,y,k,true.mean, ...){
  if (inherits(prior,"powerprior")) return(BMSEpp(prior,y=y,k=k,true.mean=true.mean, ...))
  else{
    if (ncol(prior)==1) return(BMSEuni(prior,y=y,k=k,true.mean=true.mean, ...))
    else                return(BMSEmix(prior,y=y,k=k,true.mean=true.mean, ...))
  }
}
      

############
# BMSEuni
  BMSEuni <- function(prior,...) UseMethod("BMSEuni")
  BMSEuni.default <- function(prior, ...) warning("BMSEuni does not know how to handle object")
  
  BMSEuni.betaMix <-
  function(prior,y,k,true.mean, ...){
    mean=(prior["a",1]+y)/(prior["a",1]+prior["b",1]+k)
    var=(((prior["a",1]+y)*(prior["b",1]+k-y))/((prior["a",1]+prior["b",1]+k)^2*(prior["a",1]+prior["b",1]+k+1)))
    var  + (mean-true.mean)^2
  }
  
  BMSEuni.normMix <-
  function(prior,y,k,true.mean, ...){
    sigma2=RBesT::sigma(prior)^2
    var=  1/(1/prior["s",1]^2+k/sigma2)
    mean= var*(prior["m",1]/(prior["s",1]^2)+y/(sigma2/k))
    var  + (mean-true.mean)^2
  }


############
# BMSEpp
  BMSEpp <- function(prior,...) UseMethod("BMSEpp")
  BMSEpp.default <- function(prior, ...) warning("BMSEpp does not know how to handle object")
  
  BMSEpp.normMix <-
  function(prior,y,k,true.mean, ...){
  
    sigma=RBesT::sigma(prior)
    ds=powerparameter(prior,m=y,n=k,sigma=sigma)
  
    var=  1/(1/(prior["s",1]^2/ds)+k/sigma^2)
    mean= var*(prior["m",1]/(prior["s",1]^2/ds)+y/(sigma^2/k))
    var  + (mean-true.mean)^2
  }
  
  BMSEpp.betaMix <-
  function(prior,y,k,true.mean, ...){
  #  pp=postmix.powerprior(prior,r=y,n=k, ...)
  #  pp=powerprior(prior,X=y,N=k, ...)
    pp=powerprior(prior,r=y,n=k, ...)
    BMSEuni(pp,y,k,true.mean)
  }


############
# BMSEmix
  BMSEmix <- function(prior,...) UseMethod("BMSEmix")
  BMSEmix.default <- function(prior, ...) warning("BMSEmix does not know how to handle object")
  
  
  BMSEmix.betaMix <-
  function(prior,y,k,true.mean, ...){
    posterior=postmix(prior,r=y,n=k) 
    ss=summary(posterior,probs=NULL)
    var=  ss["sd"]^2
    mean= ss["mean"]
    var  + (mean-true.mean)^2
  }
  
  BMSEmix.normMix<-
  function(prior,y,k,true.mean, ...){
    posterior= my.postmix.normMix(prior, m=y,n=k,se=RBesT::sigma(prior)/sqrt(k))
  
    sapply(posterior,function(z){
      ss= my.summary.normMix(z)
      var=  ss["sd"]^2
      mean= ss["mean"]
      var  + (mean-true.mean)^2
    })
  }


  
  
  
my.summary.normMix=function(object){
    p <- object[1, ]
    m <- object[2, ]
    v <- object[3, ]^2
    m2 <- v + m^2
    mmix <- sum(p * m)
    vmix <- sum(p * (m2 - (mmix)^2))

    c(mean = mmix, sd = sqrt(vmix))
}
# BMSE.norm.mix <-
# function(y,k,true.mean,sigma,PRIOR){
#   posterior= my.postmix.normMix(PRIOR, m=y,n=k,se=sigma/sqrt(k))
#   integrate(function(x) dmix(posterior,x)*(x-true.mean)^2,-Inf,Inf)$value
# }
# 
# BMSE.norm.mix <-
# function(y,k,true.mean,sigma,PRIOR){
#   posterior= my.postmix.normMix(PRIOR, m=y,n=k,se=sigma/sqrt(k))
# 
#   ss=summary(posterior,probs=NULL)
#   var=  ss["sd"]^2
#   mean= ss["mean"]
#   var  + (mean-true.mean)^2
# }
# BMSE.norm.mix <-
# function(y,k,true.mean,sigma,PRIOR){
#   posterior= my.postmix.normMix(PRIOR, m=y,n=k,se=sigma/sqrt(k))
# 
#   sapply(posterior,function(z){
#     ss= my.summary.normMix(z)
#     var=  ss["sd"]^2
#     mean= ss["mean"]
#     var  + (mean-true.mean)^2
#   })
# }


# BMSE.binom.mix <-
# function(y,k,true.mean,PRIOR){
#   posterior=postmix(PRIOR,r=y,n=k) #
#   integrate(function(x) dmix(posterior,x)*(x-true.mean)^2,0,1)$value
# }
# BMSE.binom.mix <-
# function(y,k,true.mean,PRIOR){
#   posterior=postmix(PRIOR,r=y,n=k) 
#   ss=summary(posterior,probs=NULL)
#   var=  ss["sd"]^2
#   mean= ss["mean"]
#   var  + (mean-true.mean)^2
# }


