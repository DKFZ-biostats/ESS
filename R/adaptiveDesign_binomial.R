adaptiveDesign_binomial <- # doesn't support powerprior because oc2S() used to compute power doesn't support them
  function(ctl.prior, treat.prior, 
           N1, Ntarget, Nmin, M, 
           pc, pt, 
           discard.prior=TRUE,vague=mixbeta(c(1, 1,1)),
           ess="ecss",
           ehss.method="mix.moment", subtractESSofVague=TRUE,
           min.ecss, D=MSE,
           decision) {
  P <- try(data.frame(pc = pc, pt = pt))
  if (inherits(P, "try-error")) {
      stop("pc and pt need to be of same size")
  }
  power_all <- matrix(0, N1+1, nrow(P))

# calculate the different possible ESS values which we get after stage1
  ESSfirst=ESSstage1=N <- vector("double", N1+1)
  ctl.prior.use=list()
  for(r in 0:N1){
      # ESS 
        if (ess=="ehss"){
          post=suppressMessages(postmix(ctl.prior, r=r/N1*Ntarget, n=Ntarget))
          ESSfirst[r+1]=ESSstage1[r+1] <- round(suppressMessages(suppressWarnings(ess(post, method=ehss.method))),0)-Ntarget
          if (subtractESSofVague) ESSfirst[r+1]=ESSstage1[r+1]=ESSstage1[r+1]-round(suppressMessages(suppressWarnings(ess(vague, method=ehss.method))),0)
        } else if (ess=="ecss"){
          U=ecss(ctl.prior,r=r,n=N1,
               n.target=Ntarget, min.ecss=min.ecss ,prior.base=vague,D=D
               )$ecss
          ESSfirst[r+1]=ESSstage1[r+1]  <- round(unname(U),0)
        } else stop("ess must be either ehss or ecss")
      
      # Discarding
        if (ESSstage1[r+1]<0 & discard.prior==TRUE){#use vague prior in case of conflict
          ctl.prior.use[[r+1]]=vague
          # ess of vague
          if (ess=="ehss"){
            post=suppressMessages(postmix(ctl.prior.use[[r+1]], r=r/N1*Ntarget, n=Ntarget))
            if (!subtractESSofVague){
              ESSstage1[r+1] <- round(suppressMessages(suppressWarnings(ess(post, method=ehss.method))),0)-Ntarget  
            } else ESSstage1[r+1] <- 0
            
          } else if (ess=="ecss"){
            U=ecss(ctl.prior.use[[r+1]],r=r,n=N1,
                 n.target=Ntarget, min.ecss=min.ecss ,prior.base=vague,D=D
                 )$ecss
            ESSstage1[r+1]  <- round(unname(U),0)
          }        
        } else ctl.prior.use[[r+1]]=ctl.prior

    # number of patients enrolled in stage 2
      N2 <- pmax(Ntarget-(ESSstage1[r+1]+N1), Nmin) 
    
    # total number of patients enrolled
      N[r+1] <- N1 + N2
  
    # calculate for each scenario and sample size of the control the power
      design_calc <- oc2S(treat.prior, ctl.prior.use[[r+1]], M, N[r+1], decision)
      power_all[r+1,] <- design_calc(P$pt, P$pc)
  }

# finally take the mean with the respective weight which corresponds
# to the weight how the respective sample size occur
  w <- sapply(P$pc, function(p) dbinom(0:N1, N1, p))

pdat= data.frame(power=colSums(power_all * w), samp=colSums(w * N),ESS=colSums(w * ESSfirst),ESSwithDiscarding=colSums(w * ESSstage1))

###################
#MSE
Nmax <- max(N)

## calculate for each i = 0 to N1 possible responders in stage one
## the posterior when observing 0 to N[i] in total. Calculate for
## each scenario outcome E(p) and E(p^2)
m  <- matrix(0, N1+1, Nmax+1)
m2 <- matrix(0, N1+1, Nmax+1)
for(r1 in 0:(N1)) {
    for(r2 in 0:N[r1+1]) {
        res <- summary(postmix(ctl.prior.use[[r1+1]],r=r2,n=N[r1+1]))[c("mean", "sd")]
        m[r1+1,r2+1] <- res["mean"]
        m2[r1+1,r2+1] <- res["sd"]^2 + m[r1+1,r2+1]^2
    }
}

## now collect the terms correctly weighted for each assumed true rate
bias <- rMSE <- c()
pt <- seq(0,1,length=101)
for(p in pt) {
    ## weight for each possible N at stage 1
    wnp <- dbinom(0:N1, N1, p)

    ## E(p) and E(p^2) for each possible N at stage 1
    Mnm <- rep(0, N1+1)
    Mnm2 <- rep(0, N1+1)

    ## for a given weight at stage 1....
    for(r1 in 0:N1) {
        ## weights of possible outcomes when having n draws in
        ## stage1, we go up to Nmax+1 to get a vector of correct
        ## length; all entries above n are set to 0 from dbinom as
        ## expected as we can never observe more counts than the
        ## number of trials...
        wp <- dbinom(0:Nmax, N[r1+1], p)

        Mnm[r1+1]  <- sum(m[r1+1,]  * wp)
        Mnm2[r1+1] <- sum(m2[r1+1,] * wp)
    }

    ## ... which we average over possible outcomes in stage 1
    Mm  <- sum(wnp * Mnm)
    Mm2 <- sum(wnp * Mnm2)

    bias <- c(bias, (Mm - p))
    rMSE <- c(rMSE, sqrt(Mm2 - 2 * p * Mm + p^2))
}
msedat= data.frame(p=pt, bias=bias, rMSE=rMSE)

list(power=pdat,MSE=msedat, prior=ctl.prior)


}
