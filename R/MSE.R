MSE <- function(prior,...) UseMethod("MSE")
MSE.default <- function(prior, ...) warning("MSE does not know how to handle object")


MSE.mix <-function(prior,y,k,true.mean, ...){
  if (inherits(prior,"powerprior")) return(MSEpp(prior,y=y,k=k,true.mean=true.mean, ...))
  if (inherits(prior,"commensurateEB")) return(MSEcommensurateEB(prior,y=y,k=k,true.mean=true.mean, ...))
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
# MSEcommensurate
  MSEcommensurateEB <- function(prior,...) UseMethod("MSEcommensurateEB")
  MSEcommensurateEB.default <- function(prior, ...) warning("MSEcommensurateEB does not know how to handle object")
  
MSEcommensurateEB.normMix <-
  function(prior,y,k,true.mean,l.eb=0,u.eb=Inf, ...){
  
    sigma=RBesT::sigma(prior)
  
    # MMLE of tau, Hobbs et al 2012 (5) ##
      Delta.hat <- y -  prior["m",]
      tau.EB <- 1/max( min( (Delta.hat^2)-(sigma^2/k)-prior["s",1]^2, u.eb ), l.eb )
    var=  1/(1/(prior["s",1]^2+(1/tau.EB)+ k/sigma^2) ) 
    mean=var*(prior["m",1]/(prior["s",1]^2+(1/tau.EB)) + y/(sigma^2/k)   )
    (mean-true.mean)^2
   
}
MSEcommensurateEB.betaMix <-
  function(prior,y,k,true.mean, ...){
  
    r=y
    n=k
    n0=prior["a",]+prior["b",]
    y0=prior["a",]
      
      k1=seq(0,500,by=.2)
      marg=rep(NA,length(k1))
      
      for (i in 1:length(k1)) {
        integrandU.V=Vectorize(function(th0C){
          integrand=function(th,th0=th0C,k=k1[i]) dbinom(r,n,th)*dbinom(y0,n0,th0)*dbeta(th,k*th0,k*(1-th0))
          as.numeric(integrate(integrand,0.001,0.999)[1])
        })
        marg[i]=as.numeric(integrate(integrandU.V,0.001,0.999)[1])
      }
      
      k_emp=k1[which.max(marg)]
      marg_e=marg[which.max(marg)]
      
      integrandU.V=Vectorize(function(thN){
        integrand=function(th0,th=thN,k=k_emp)  dbinom(r,n,th)*dbinom(y0,n0,th0)*dbeta(th,k*th0,k*(1-th0))
        as.numeric(integrate(integrand,0,1)[1])/marg_e
      })
      
      xg=seq(0.001,0.999,by=0.001)
      iU=integrandU.V(xg)
      mean.thC=sum(xg*iU/sum(iU))
   #   var.thC=sum(xg^2*iU/sum(iU))-mean.thC^2

      (mean.thC-true.mean)^2
  
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
