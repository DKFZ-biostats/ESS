plot.ecssReimherr <-
function(x,n.target,type="target",power=NULL,tol=1e-5,...){
  ESS=ecssReimherr(x,n.target,type=type,power=power,tol=tol)
  Ub=splinefun(x=seq_along(x$U.base), y=x$U.base)
  Ui=splinefun(x=seq_along(x$U.info), y=x$U.info)
  
  plot(x$U.base,type="l",col=1,lwd=2,lty="dashed",sub="root of Uinfo(k)-Ubase(k+ESS)",...)
  lines(x$U.info,type="l",col=2)
  # curve(fb,lty="dashed",add = T,lwd=2)
  # curve(ft,lty="dashed",add = T,col=2,lwd=2)
  #abline(h=Ui(k),col=2)
  arrows(x0=n.target, y0=Ui(n.target), x1 = n.target+ESS, y1 = Ui(n.target),length=.15)
  text(n.target+ESS/2,Ui(n.target),paste0("ESS=",round(ESS,1)),cex=1,font=3,pos=3)
  
  abline(v=n.target,col=2)
  abline(v=n.target+ESS,col=1,lwd=2,lty="dashed")
  legend("topright",legend = c("base","info"),text.col = c(1,2))
  return(c(ESS=ESS))
}
