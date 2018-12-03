adaptiveDesign_normal <-
function(ctl.prior, treat.prior,
         N1, Ntarget, Nmin, M, 
         muc, mut,sc,st, sc.known=TRUE,st.known=TRUE,
         discard.prior=TRUE, vague=mixnorm(c(1, 0,10)),
         ess="ecss",ehss.method="mix.moment", min.ecss, D=MSE,
         decision, 
         nsim=100,cores=1, seed=123, progress="text"
         ) {

  seedL=sapply(gather(nsim,seed=seed),function(x) x[2:7])

  
  MU <- try(data.frame(muc = muc, mut = mut))
  if (inherits(MU, "try-error")) {
      stop("muc and mut need to be of same size")
  }
    #   cl=makeForkCluster(cores)
    # registerDoParallel(cl, cores=cores)

    # res=foreach(z=1:nrow(MU))%:%
    #                foreach(it=1:nsim) %dopar%  simDesign(seed=seedL[,it],
                                                           # ctl.prior=ctl.prior, treat.prior=treat.prior, 
                                                           # N1=N1, Ntarget=Ntarget, Nmin=Nmin, M=M,
                                                           # muc=MU$muc[z],mut=MU$mut[z],sc=sc,st=st,
                                                           # discard.prior=discard.prior,
                                                           # ess=ess,ehss.method=ehss.method, min.ecss=min.ecss,
                                                           # decision=decision)
    # stopCluster(cl)
  requireNamespace("plyr", quietly = TRUE)
  res=plyr::llply(1:nrow(MU), function(z)
                    mclapply(1:nsim,function(it) simDesign_normal(seed=seedL[,it],
                                                           ctl.prior=ctl.prior, treat.prior=treat.prior,
                                                           N1=N1, Ntarget=Ntarget, Nmin=Nmin, M=M,
                                                           muc=MU$muc[z],mut=MU$mut[z],sc=sc,st=st, sc.known=sc.known,st.known=st.known,
                                                           discard.prior=discard.prior, vague=vague,
                                                           ess=ess,ehss.method=ehss.method, min.ecss=min.ecss, D=D,
                                                           decision=decision),
                             mc.cores = cores), .progress=progress
                                      )
  # res=mclapply(1:nrow(MU), function(z)
  #                   lapply(1:nsim,function(it) simDesign_normal(seed=seedL[,it],
  #                                                          ctl.prior=ctl.prior, treat.prior=treat.prior, 
  #                                                          N1=N1, Ntarget=Ntarget, Nmin=Nmin, M=M,
  #                                                          muc=MU$muc[z],mut=MU$mut[z],sc=sc,st=st, sc.known=sc.known,st.known=st.known,
  #                                                          discard.prior=discard.prior, vague=vague,
  #                                                          ess=ess,ehss.method=ehss.method, min.ecss=min.ecss, D=D,
  #                                                          decision=decision)),
  #                            mc.cores = cores
  #                                     )

   res2=lapply(1:length(res),function(z) rbind(muc=MU$muc[z],mut=MU$mut[z],simplify2array(res[[z]])))
   names(res2)=apply(MU,1,paste,collapse="vs.")

  list(list=res2)
}


simDesign_normal <-
function(seed,
         ctl.prior, treat.prior,
         N1, Ntarget, Nmin, M,
         muc,mut,sc,st, sc.known=TRUE,st.known=TRUE,
         discard.prior=TRUE,vague=mixnorm(vague=c(1, 0,10)),
         ess="ecss",ehss.method="mix.moment", min.ecss, D=BMSE,
         decision
        ){
  set.seed(seed,"L'Ecuyer-CMRG")
  # interim
    obs.interim=rnorm(N1,muc,sc)
    if (is.logical(sc.known)){
      if (sc.known) sigma.c=sc else sigma.c=sd(obs.interim)
    } else if (is.numeric(sc.known)) sigma.c=sc.known
   
    
  # ESS 
    if (ess=="ehss"){
      post=suppressMessages(postmix(ctl.prior, n=Ntarget,m=mean(obs.interim),se=sigma.c/sqrt(Ntarget)))
      ESSstage1=ESSfirst<- round(suppressMessages(suppressWarnings(ehss(post, method=ehss.method))),0)-Ntarget
    } else if (ess=="ecss"){
      U=ecss(ctl.prior,m=mean(obs.interim),sigma=sigma.c, n=N1,
           n.target=Ntarget, min.ecss=min.ecss ,prior.base=vague,D=D
           )$ecss
      ESSstage1=ESSfirst <- round(unname(U),0)
    } else stop("ess must be either ehss or ecss")
  
  # Discarding
    if (ESSstage1<0 & discard.prior==TRUE){#use vague prior in case of conflict
      ctl.prior.use=vague
      ESSstage1=0
    } else ctl.prior.use=ctl.prior

  # number of patients enrolled in stage 2
    N2 <- pmax(Ntarget-(ESSstage1+N1), Nmin) 
  # total number of patients enrolled
    N <- N1 + N2

  # final analysis
    ctl.obs.final=c(obs.interim,rnorm(N2,muc,sc))
    treat.obs.final=rnorm(M,mut,st)
    
    
    if (is.logical(st.known)){
      if (st.known) sigma.t=st else sigma.t=sd(treat.obs.final)
    } else if (is.numeric(st.known)) sigma.t=st.known
    post.t=suppressMessages(postmix(treat.prior, n=M,m=mean(treat.obs.final),se=sigma.t/sqrt(M)))

    if (is.logical(sc.known) & sc.known==FALSE){
      post.c=suppressMessages(postmix(ctl.prior.use,data=ctl.obs.final))
    } else {
      post.c=suppressMessages(postmix(ctl.prior.use, n=N,m=mean(ctl.obs.final),se=sigma.c/sqrt(N)))
    }

 
    rej=decision( post.t,  post.c) 
    prob=pmixdiff(post.t, post.c, 0, lower.tail = FALSE)
    mse=  integrate(function(x) dmixdiff(post.t,post.c,x)*(x-(mut-muc))^2,-Inf,Inf)$value #mse=  mean((rmixdiff(post.t, post.c,1e4)  -(mut-muc))^2)
    muc.post=summary(post.c,probs=NULL)["mean"]
    mut.post=summary(post.t,probs=NULL)["mean"]
    bias= (mut.post-muc.post ) -(mut-muc)
    biasSq=bias^2    

  c(power=rej,N=N,ESS=ESSfirst,mse=mse,bias=bias,biasSq=biasSq,prob=prob,
    muc.obs.interim=mean(obs.interim),
    muc.obs.final=mean(ctl.obs.final),
    mut.obs.final=mean(treat.obs.final),
    muc.post.final=muc.post,
    mut.post.final=mut.post
  )
}
