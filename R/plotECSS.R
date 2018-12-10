plotECSS=function(priorlist, grid, n.target, min.ecss,...){
 
  ECSS=lapply(priorlist, function(prior) sapply(grid, function(j) ecss(prior,true.mean=j,
                                               min.ecss=min.ecss,n.target=n.target,...)$ecss)) 
  dat=data.frame(grid=rep(grid,length(priorlist)),
                 ECSS=unlist(ECSS),
                 prior=ordered(rep(names(priorlist),each=length(grid)),levels=names(priorlist) )
                 )

  p=ggplot()+geom_line(aes(grid, ECSS, colour=prior, linetype=prior), data=dat)
  res=list(plot=p,data=dat)
  class(res)="ECSSplot"
  res
}

print.ECSSplot<-function(x, ...) print(x$plot)
