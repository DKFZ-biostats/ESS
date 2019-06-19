as.commensurateEB=function(prior){
  if (ncol(prior)!=1) stop("mixture object may only have one component")
 # if (inherits(prior,"betaMix")) attr(prior,"p.prior")=c(p.prior.a,p.prior.b)
  class(prior)=c("commensurateEB",class(prior))
  prior
}



################
# postmix
  postmix.commensurateEB=function(priormix, data, n, r, m, se, l.eb=0,u.eb=Inf, ...){
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

      
      Delta.hat <- m -  priormix["m",]
      ## MMLE of tau, Hobbs et al 2012 (5) ##
#      v0=priormix["s",]^2
#        tau.EB <- 1/max( min( (Delta.hat^2)-se^2-v0, u.eb ), l.eb )
#      s.post= sqrt( 1/( 1/se^2 + 1/(v0+(1/tau.EB)) ) ) 
        tau.EB.inv <- max( min( (Delta.hat^2)-se^2-priormix["s",]^2, u.eb ), l.eb )
      s.post= sqrt( 1/( 1/se^2 + 1/(priormix["s",]^2+tau.EB.inv) ) ) 
      # m.post=( m/se^2  + priormix["m",]/(v0+(1/tau.EB)) )/
      #             ( 1/se + (1/(v0+(1/tau.EB))) )
      m.post=( m/se^2  + priormix["m",]/(priormix["s",]^2+tau.EB.inv) )*s.post^2
      return(mixnorm(c(1,m.post,s.post), sigma=se * sqrt(n)))
    }
    if (inherits(priormix,"betaMix")){
      
      
      if (!missing(data)) {
          assert_that(all(data %in% c(0, 1)))
          r <- sum(data)
          n <- length(data)
      }
      # if (missing(p.prior.a)| missing(p.prior.b)){
      #   if (!is.null(attr(priormix,"p.prior"))){
      #     p.prior.a=attr(priormix,"p.prior")[1]
      #     p.prior.b=attr(priormix,"p.prior")[2]
      #   } else{
      #     p.prior.a=p.prior.b=1
      #   }
      # }
      # 
      
      n0=priormix["a",]+priormix["b",]
      y0=priormix["a",]
      
#      k1=seq(0,1000,by=1)
      k1=seq(0,500,by=.1)
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
      var.thC=sum(xg^2*iU/sum(iU))-mean.thC^2
      
      # appoximates distribution by beta distribution, seems to work well, but needs to be further checked
      estBetaParams <- function(mu, var) {
        alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
        beta <- alpha * (1 / mu - 1)
        return(params = list(alpha = alpha, beta = beta))
      }
      ab=estBetaParams(mean.thC,var.thC)
      mixbeta(c(1,ab$alpha,ab$beta))

      
      
    }
  }

# nu <- function(M0,a,b.tau,Bd){
#  return( max( min( ( (M0/a)^2 ) - b.tau, 1/Bd[1] ), 1/Bd[2] ) )
# }
# 
# ## Posterior mean for lambda|tau ## (10)
# M <- function(tau,A,B,n,n1,sig2,v0){ return(  ( A/( (n/n1)*((sig2/n)+v0+(1/tau)) ) + ( B/(sig2/n1)) ) /
#                                                 ( 1/V(tau,n,n1,sig2,v0) ) 
#                                             ) }
# ## Posterior Variance of lambda|tau ## (9)
# V <- function(tau,n,n1,sig2,v0){ return( 1/( 1/(((n/n1)^2)*((sig2/n)+v0+(1/tau))) + ( ( n1*(1-(n1/n)) )/sig2 ) ) ) }
# 
# 
# gen.Data <- function(pars){
#  data. <- list(A=rnorm(1, pars$Ex.A, sqrt(pars$V.A)),
#                B=rnorm(1, pars$Ex.B, sqrt(pars$V.B)) )
#  data.$M0 <- data.$A*pars$a - data.$B/( 1-(1/pars$a) )
#  return(data.)
# }
# 
# eb <- function(pars,Bd){
#  data. <- gen.Data(pars)
#  tau <- 1/nu(data.$M0,pars$a,pars$b.tau,Bd)
#  EHSS <- pars$sig2/(pars$v0+(1/tau))
#  return( c(M(tau,data.$A,data.$B,pars$n,pars$n1,pars$sig2,pars$v0), ## Post Mean lambda
#            V(tau,pars$n,pars$n1,pars$sig2,pars$v0), ## Post var lambda
#            EHSS ## EHSS
#            ) )
# }
