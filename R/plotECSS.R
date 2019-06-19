plotECSS=function(priorlist, grid, n.target, min.ecss,progress = "text",...){
#  utils::globalVariables(c("prior"))
  ECSS=suppressMessages(llply(priorlist, function(prior) sapply(grid, function(j) ecss(prior,true.mean=j,
                                               min.ecss=min.ecss,n.target=n.target,...)$ecss),.progress = progress) )
  dat=data.frame(grid=rep(grid,length(priorlist)),
                 ECSS=unlist(ECSS),
                 prior=ordered(rep(names(priorlist),each=length(grid)),levels=names(priorlist) )
                 )

  p=ggplot()+geom_line(aes_string("grid", "ECSS", colour="prior"), data=dat)+
             geom_line(aes(grid, ECSS), colour="white",linetype="dashed", data=subset(dat,dat$ECSS==min.ecss )) +
             geom_line(aes(grid, ECSS), colour="white",linetype="dashed", data=subset(dat,dat$ECSS==n.target)) 
  res=list(plot=p,data=dat)
  class(res)="ECSSplot"
  res
}

print.ECSSplot<-function(x, ...) print(x$plot)
