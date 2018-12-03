plot.ecssReimherr.binom <-
function(obj,n.target,...){
  ESS=ecssReimherr.binom(obj,n.target)
    Ub=approxfun(x=seq_along(obj$U.base), y=obj$U.base,method="const",f=0)
    Ui=approxfun(x=seq_along(obj$U.info), y=obj$U.info,method="const",f=0)
   
     plot(obj$U.base,type="l",col=1,lwd=2,lty="dashed",sub="root of Uinfo(n.target)-Ubase(n.target+ESS)",...)
  lines(obj$U.info,type="l",col=2)
   # curve(Ub,lty="dashed",add = T,lwd=2)
   # curve(Ui,lty="dashed",add = T,col=2,lwd=2)
  abline(h=Ui(n.target),col=2)
  abline(v=n.target,col=2)
  abline(v=n.target+ESS,col=1,lwd=2,lty="dashed")
  legend("topright",legend = c("base","info"),text.col = c(1,2))

   
  ESS
}
