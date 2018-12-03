MSE <- function(prior,...) UseMethod("MSE")
MSE.default <- function(prior, ...) warning("MSE does not know how to handle object")


MSE.mix <-function(prior,y,k,true.mean, ...){
  if (inherits(prior,"powerprior")) return(MSEpp(prior,y=y,k=k,true.mean=true.mean, ...))
  else{
    if (ncol(prior)==1) return(MSEuni(prior,y=y,k=k,true.mean=true.mean, ...))
    else                return(MSEmix(prior,y=y,k=k,true.mean=true.mean, ...))
  }
}
      

############
# MSEuni
MSEuni <- function(prior,...) UseMethod("MSEuni")
MSEuni.default <- function(prior, ...) warning("MSEuni does not know how to handle object")

MSEuni.betaMix <-
  function(prior,y,k,true.mean, ...){
    mean=(prior["a",1]+y)/(prior["a",1]+prior["b",1]+k)
 #   var=(((prior["a",1]+y)*(prior["b",1]+k-y))/((prior["a",1]+prior["b",1]+k)^2*(prior["a",1]+prior["b",1]+k+1)))
    (mean-true.mean)^2
  }

MSEuni.normMix <-
function(prior,y,k,true.mean, ...){
  sigma2=RBesT::sigma(prior)^2
  var=  1/(1/prior["s",1]^2+k/sigma2)
  mean= var*(prior["m",1]/(prior["s",1]^2)+y/(sigma2/k))
  (mean-true.mean)^2
}


############
# MSEpp
  MSEpp <- function(prior,...) UseMethod("MSEpp")
  MSEpp.default <- function(prior, ...) warning("MSEpp does not know how to handle object")
  
  MSEpp.normMix <-
  function(prior,y,k,true.mean, ...){
  
    sigma=RBesT::sigma(prior)
    ds=powerparameter(prior,m=y,n=k,sigma=sigma)
  
    var=  1/(1/(prior["s",1]^2/ds)+k/sigma^2)
    mean= var*(prior["m",1]/(prior["s",1]^2/ds)+y/(sigma^2/k))
    (mean-true.mean)^2
 }
  
  MSEpp.betaMix <-
  function(prior,y,k,true.mean, ...){
  #  pp=postmix.powerprior(prior,r=y,n=k, ...)
  #  pp=powerprior(prior,X=y,N=k, ...)
    pp=powerprior(prior,r=y,n=k, ...)
    MSEuni(pp,y,k,true.mean)
  }

############
# MSEmix
  MSEmix <- function(prior, ...) UseMethod("MSEmix")
  MSEmix.default <- function(prior, ...) warning("MSEmix does not know how to handle object")

  MSEmix.betaMix <-
  function(prior,y,k,true.mean, ...){
    posterior=postmix(prior,r=y,n=k) 
    ss=summary(posterior,probs=NULL)
    #var=  ss["sd"]^2
    mean= ss["mean"]
    (mean-true.mean)^2
  }
  
  MSEmix.normMix<-
  function(prior,y,k,true.mean, ...){
    posterior= my.postmix.normMix(prior, m=y,n=k,se=RBesT::sigma(prior)/sqrt(k))
  
    sapply(posterior,function(z){
      p <- z[1, ]
      m <- z[2, ]
      mean= sum(p * m)
      (mean-true.mean)^2
    })
  }


# MSEmix.normMix <-
# function(prior,y,k,true.mean, ...){
#   posterior= my.postmix.normMix(prior, m=y,n=k,se=RBesT::sigma(prior)/sqrt(k))
#   mean=summary(posterior,probs=NULL)["mean"]
#   (mean-true.mean)^2
# }
