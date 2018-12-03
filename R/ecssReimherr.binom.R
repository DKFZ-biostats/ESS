ecssReimherr.binom= Vectorize(function(obj,n.target,power=NULL,tol=1e-3) {
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
  
   if(abs(obj$U.base[n.target]-obj$U.info[min(n.target,length(obj$U.info))])<tol |
      abs(obj$U.base[n.target]-obj$U.info[n.target])>(1-tol)
      ){#curves are practically overlapping (maybe we need to add the case when one is always 0 and the other always 1?)
   return(c(ESS=0))
  }   else {
      if ((baseBelowInfo &                             #no Ess solution
           baseIncreasing==FALSE)                          #Base is decreasing (type I/MSE)
          | (baseAboveInfo &                          #no Ess solution   
             baseIncreasing==TRUE)                     #Base is increasing (Power)
         ){
          
      return(c(ESS=-n.target))}else {
        if ( (baseBelowInfo
             & baseIncreasing==TRUE)              #Base is increasing (Power)
            | (baseAboveInfo 
              & baseIncreasing==FALSE )           #Base is decreasing (type I/MSE)
            ){
          
        return(c(ESS=ReimMax-n.target))
          } else {
    Ub=approxfun(x=seq_along(obj$U.base), y=obj$U.base,method="const",f=0)
    #Ui=approxfun(x=seq_along(obj$U.info), y=obj$U.info,method="const",f=0)

#   n.target-optimize(function(xx) abs(fb(x)-ft(xx)) , interval=c(1,length(obj$U.base)) )$minimum)
      # ESS=min(optimize(function(xx) abs(Ui(n.target)-Ub(n.target+xx)) , interval=c(-n.target+1,length(obj$U.base)-n.target) )$minimum)
      # if (abs(Ui(n.target)-Ub(n.target+ESS))>tol) ESS=NA
  #    xx=seq(-n.target+1,length(obj$U.base)-n.target,length.out=1e4)
      xx=seq(-n.target+1,length(obj$U.base)-n.target,length.out=length(obj$U.base))
      ESS=xx[which.min(abs(obj$U.info[n.target]-Ub(n.target+xx)))]
       return(c(ESS=ESS))
        }
      }
    }
    },vectorize.args="n.target")
