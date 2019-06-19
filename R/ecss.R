ecss <- function(prior,...) UseMethod("ecss")
ecss.default <- function(prior, ...) warning("ecss does not know how to handle object")
print.ecss<-function(x, ...) print(c(ecss=x$ecss), ...)


ecss.mixture.prior <- function(prior,...){
  prior=as.mix(prior)
  ecss(prior,...)
}

ecss.mix <-
function(prior,
         data,
         n, r,
         m, sigma, se,
         true.mean,
         n.target, min.ecss,# vector of n.target allowed 
         prior.base,  
         D=MSE,
         by,
         grid.length=50, min.q,
         cores=1,
         integrate=FALSE, subdivisions=100L,...){
  
  call=match.call(expand.dots = T)
  n.max=max(n.target)-min.ecss
args=as.list(match.call(expand.dots = FALSE))
  
  if (inherits(prior,"normMix") | inherits(prior,"normFix")){
    if (missing(prior.base)){
      warning("Setting prior.base to N(0,100)")
      prior.base=mixnorm(vague=c(1, 0,10), sigma=sigma)
    }
    if (max(prior["s",])>prior.base["s",]) warning("At least one mixture component is more dispersed than prior.base. Consider increasing dispersion of prior.base.")
    if (!inherits(prior.base,"normMix") ) stop("Object prior.base must be a normMix object.")
    if (ncol(prior.base)!=1) stop("prior.base may only have one mixture component")

    if (missing(by)) by=5
   if (missing(min.q)) min.q=1e-6

    if (!missing(true.mean)){
      data.mean=true.mean  
    } else {
      if (!missing(data)) {
          m <- mean(data)
          sigma=sd(data)
          n=length(data)
      }
      else {
          if (missing(m) & (missing(sigma) &((missing(n) | missing(se))) ))
              stop("Either raw data or summary data (m and se) must be given.")
          if (!missing(se) & !missing(n) &  missing(sigma)) {
              sigma <- se * sqrt(n)
              # message(paste0("Updating reference scale to ", sigma(priormix), 
              #     ".\nIt is recommended to use the sigma command instead.\nSee ?sigma or ?mixnorm."))
          }
          if (!missing(se) & missing(n) &  !is.null(RBesT::sigma(prior))) {
              message("Using default prior reference scale ", RBesT::sigma(prior))
              sigma=RBesT::sigma(prior)
              n=(sigma/se)^2
          }

          if (missing(se) & missing(sigma) & !missing(n)&  !is.null(RBesT::sigma(prior))) {
              message("Using default prior reference scale ", RBesT::sigma(prior))
              sigma <- RBesT::sigma(prior)
          }
       }
  
  #    data.mean=m
      
      sigma2.post=(1/prior.base["s",1]^2+n/sigma^2)^(-1)
      data.mean= sigma2.post*(prior.base["m",1]/(prior.base["s",1]^2)+m/(sigma^2/n))
    }
    
    sigma(prior)=sigma(prior.base)=sigma
    #sigma2=RBesT::sigma(prior)^2
 
    weight=function(y,k) dnorm(y,data.mean,sigma/sqrt(k))
#    weight=function(y,k) dnorm(y,post.mean,sqrt(sigma^2+sigma2.post)/sqrt(k))
  #  weight=function(y,k) dnorm(y,post.mean,sqrt(sigma^2/k+sigma2.post))
  } else if (inherits(prior,"betaMix")){
    if (missing(prior.base)){
      warning("Setting prior.base to uniform")
      prior.base=mixbeta(               c(1.0, 1, 1))
    }
    if (!inherits(prior.base,"betaMix") ) stop("Object prior.base must be a betaMix object.")
    if (ncol(prior.base)!=1) stop("prior.base may only have one mixture component")
    
    if (inherits(prior,"powerprior")){
      if (!all(prior.base[c("a","b"),1]==attr(prior,"p.prior"))) stop("Parameters p.prior.a and p.prior.b in as.powerprior() must be equal to parameters in prior.base.")
    }
    
    if (missing(by)) by=1
   if (missing(min.q)) min.q=0
    
    if (!missing(true.mean)) data.mean=true.mean else {
      if (!missing(data)) {
        assert_that(all(data %in% c(0, 1)))
        r <- sum(data)
        n <- length(data)
      }
#     data.mean=r/n 
      shape1.post=(r+prior.base["a",1])
      shape2.post=n-r+prior.base["b",1]
      data.mean= shape1.post/(shape1.post+shape2.post)
    }
    

   weight=function(y,k) dbinom(y,size=k,data.mean)    
 #    weight=function(y,k) dbetabinom.ab(y,size=k,shape1.post,shape2.post)    
  } else stop("prior not supported")
  
  if (!is.function(D)) stop ("D must be a function.")
  Dfun.vec=function(y,k,PRIOR) weight(y,k)*do.call(D,c(list(prior=PRIOR,y=y,k=k,true.mean=data.mean),args$...))
  Dfun=Vectorize(Dfun.vec,vectorize.args = "y")
    
  nvec=sort(union(seq(1,n.max,by=by),c(n.target,n.max)))
  res=mclapply(nvec, function(k){
    if (inherits(prior,"normMix")| inherits(prior,"normFix")){
      # U.info
        if (integrate ) U.info= integrate(Dfun,-Inf,Inf, k=k, PRIOR=prior,subdivisions = subdivisions,stop.on.error = F)$value
        else {
          nq=qnorm(min.q)*sigma/sqrt(k)
          grid=seq(data.mean+nq,data.mean-nq,length.out=grid.length)
          d.grid=diff(grid[1:2])
          U.info=sum(d.grid*Dfun.vec(grid,k=k,PRIOR=prior))
        }  
      # U.base
        if (k==1 | k%in%n.target | k==n.max){
          if (integrate) U.base= integrate(Dfun,-Inf,Inf, k=k, PRIOR=prior.base,subdivisions = subdivisions,stop.on.error = F)$value
          else {
            U.base=sum(d.grid*Dfun.vec(grid,k=k,PRIOR=prior.base))  
          } 
        } else U.base=NA
    } else if (inherits(prior,"betaMix")){
      grid=qbinom(min.q,k,data.mean):qbinom(1-min.q,k,data.mean) #0:k if min.q==0
      U.info= sum(Dfun(grid,k=k,PRIOR=prior))
      U.base=ifelse((k==1 | k%in%n.target | k==n.max),
                    yes = sum(Dfun(grid,k=k,PRIOR=prior.base)), no=NA)

          }
    list(U.base=U.base,U.info=U.info)
  },mc.cores = cores)
  U.base=sapply(res,function(z) z$U.base)
  U.info=sapply(res,function(z) z$U.info)

  res=list(call=call,U.base=U.base,U.info=U.info,n=nvec,n.target=n.target) 
  class(res)=c("ecss",paste0("ecss.",class(prior)[1]))
  ESS=sapply(n.target, function(k){
    ess=ecss(res,k)
    ess$ecss
  } )
  res=c(res,ecss=list(ESS))
  class(res)=c("ecss",paste0("ecss.",class(prior)[1]))
  res
}




ecss.ecss <-
function(prior, n.target,# type="target", 
         tol=1e-6, tol2=1e-3, ...) {
  if (!(n.target%in%prior$n)) stop("'n.target' was not evaluated in ecss")
  n.max=max(prior$n)

  infoBelowBase= (max(prior$U.info,na.rm=TRUE)<prior$U.base[which(prior$n==n.target)])
  baseBelowInfo= (prior$U.base[which(prior$n==n.target)]<min(prior$U.info,na.rm=TRUE))

  baseIncreasing= FALSE #(prior$U.base[1]<prior$U.base[length(prior$U.base)])   
  
  
  if(abs(prior$U.info[which(prior$n==n.target)]-prior$U.base[which(prior$n==n.target)])<tol) {#curves are practically overlapping (maybe we need to add the case when one is always 0 and the other always 1?)
      prior$ecss=NULL
      res=c(prior,ecss=0)
      class(res)=c("ecss")
      return(res)
  } else {
#      if ((infoBelowBase & baseIncreasing==FALSE) | (baseBelowInfo & baseIncreasing) ){
      if ((infoBelowBase & baseIncreasing==FALSE) ){
            message("Maximum estimable ECSS attained (=n.target). ECSS could likely be underestimated.")
            prior$ecss=NULL
            res=c(prior,ecss=n.target)
            class(res)=c("ecss")
            return(res)
      } else {
#        if ((infoBelowBase &   baseIncreasing) |  (baseBelowInfo & baseIncreasing==FALSE) ){
        if ((baseBelowInfo & baseIncreasing==FALSE) ){
          message("Minimum ECSS attained (=min.ecss). ECSS could be overestimated. Decrease argument min.ecss to lower negative value.")
          prior$ecss=NULL
          res=c(prior,ecss=-(n.max-n.target))
          class(res)=c("ecss")
          return(res)
        } else {
#          if (type=="target") { #evaluation at Ntarget
            Ui=splinefun(x=prior$n, y=prior$U.info)
            ESS <- min(optimize(function(xx) abs(prior$U.base[which(prior$n==n.target)]-Ui(n.target+xx)) , interval=c(-n.target+1,n.max-n.target) )$minimum)
            if (abs(prior$U.base[which(prior$n==n.target)]-Ui(n.target+ESS))>tol2) ESS=NA # might be problematic and setting ESS=NA might be better
            prior$ecss=NULL
            res=c(prior,ecss=-ESS)
            class(res)=c("ecss")
            return(res)
 #         }
          # if (type=="regression") {
          #   if (sum(!is.na(prior$U.base))==3) stop("if type='regression' you need to specify n.target=NULL")
          #   Ui=splinefun(x=prior$n, y=prior$U.info)
          #     ESS0=ecss(prior,n.target)$ecss
          #     ESSc <- sapply(1:as.integer(sqrt(n.target)), function(kk) ecss(prior,kk)$ecss)
          #     x=c(1:as.integer(sqrt(n.target)))
          #     lmEss=lm(ESSc~x)
          #   ESS<-list(c(ecss=ESS0,ecss.lm=unname(lmEss$coeff[1]+lmEss$coeff[2]*n.target)))
          #   prior$ecss=NULL
          #   res=c(prior,ecss=ESS)
          #   class(res)=c("ecss")
          #   return(res)
          # }
        }
      }
    }
}
