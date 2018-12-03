# ecssReimherr <- function(obj,...) UseMethod("ecssReimherr")
# ecssReimherr.default <- function(obj, ...) warning("ecssReimherr does not know how to handle object")
# print.ecssReimherr<-function(x, ...) print(c(ecssReimherr=x$ecss), ...)
# 

ecssReimherr= Vectorize(function(obj,n.target,type="target",power=NULL,tol=1e-5,tol2=1e-3,k0=1,K=NULL) {
  if (!(n.target%in%obj$n)) stop("'n.target' was not evaluated in ecss")
  ReimMax=length(obj$U.base)

  # baseBelowInfo=(all(obj$U.base[c(1,n.target,ReimMax)]<obj$U.info[c(1,n.target,ReimMax)]) & #base always below info and no Ess solution
  #        (max(obj$U.base,na.rm=TRUE)<(obj$U.info[n.target])))
  # baseAboveInfo= (all(obj$U.base[c(1,n.target,ReimMax)]>obj$U.info[c(1,n.target,ReimMax)]) &                      #base always above info
  #          (obj$U.info[n.target]<min(obj$U.base,na.rm=TRUE)))
  baseBelowInfo=(#all(obj$U.base[c(1,n.target,ReimMax)]<obj$U.info[c(1,n.target,ReimMax)]) & #base always below info and no Ess solution
         (max(obj$U.base,na.rm=TRUE)<(obj$U.info[n.target])))
  baseAboveInfo= (#all(obj$U.base[c(1,n.target,ReimMax)]>obj$U.info[c(1,n.target,ReimMax)]) &                      #base always above info
           (obj$U.info[n.target]<min(obj$U.base,na.rm=TRUE)))
  
  if (is.null(power)) baseIncreasing= (obj$U.base[1]<obj$U.base[ReimMax])   # TRUE for base increasing, FALSE decreasing. Can be set externally through argument "power"
  else baseIncreasing=power
  
  
#  if(sum(abs(obj$U.base-obj$U.info))<(ReimMax*tol))#curves are practically overlapping (maybe we need to add the case when one is always 0 and the other always 1?)
    if(abs(obj$U.base[n.target]-obj$U.info[n.target])<tol)#curves are practically overlapping (maybe we need to add the case when one is always 0 and the other always 1?)
      { return(c(ESS=0))}
    else {
      if ((baseBelowInfo &                             #base always below info
           baseIncreasing==FALSE)                      #Base is decreasing (type I/MSE)
          |                        
           (baseAboveInfo &               #base always above info
            baseIncreasing)                     #Base is increasing (Power)
          )
        {return(c(ESS=-n.target))}
      else {
        if ((baseBelowInfo &                    #base always below info
             baseIncreasing)                   #Base is increasing (Power)
            |   
             (baseAboveInfo &              #base always above info
              baseIncreasing==FALSE)            #Base is decreasing (type I/MSE)
            )
          {return(c(ESS=ReimMax-n.target))}
        else {
          if (type=="target") {
            #evaluation at Ntarget
            Ub=splinefun(x=seq_along(obj$U.base), y=obj$U.base)
#            Ui=splinefun(x=seq_along(obj$U.info), y=obj$U.info)
            ESS <- min(optimize(function(xx) abs(obj$U.info[n.target]-Ub(n.target+xx)) , interval=c(-n.target+1,ReimMax-n.target) )$minimum)
            #   n.target-optimize(function(xx) abs(fb(x)-ft(xx)) , interval=c(1,length(obj$U.base)) )$minimum) 
            if (abs(obj$U.info[n.target]-Ub(n.target+ESS))>tol2) ESS=NA # might be problematic and setting ESS=NA might be better
 #           attr(ESS,"diff")=obj$U.info[n.target]-Ub(n.target+ESS)
             return(c(ESS=ESS))
          }
          if (type=="regression") {
            if (sum(!is.na(obj$U.info))==3) stop("Need to specify n.target=NULL in ecss()")
            Ub=splinefun(x=seq_along(obj$U.base), y=obj$U.base)
            ESSc <- sapply(k0:K, function(kk) min(optimize(function(xx) abs(obj$U.info[kk]-Ub(kk+xx)) , interval=c(-kk+1,ReimMax-kk) )$minimum))
            lmEss=lm(ESSc~I(k0:K))
            ESS<-(sum(coef(lmEss)*c(1,n.target)))
            return(c(ESS=ESS))
          }
        }
      }
    }
},vectorize.args="n.target")

