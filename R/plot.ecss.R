plot.ecss <- function(x,n.target,#type="target",
                      tol=1e-6,tol2=1e-3,...){
    ESS=ecss(x,n.target,#type=type,
             tol=tol,tol2=tol2)$ecss

    Ub=splinefun(x=x$n, y=x$U.base)
    Ui=splinefun(x=x$n, y=x$U.info)
    
    plot(x$n,x$U.info,type="l",col=2,ylab=x$D,xlab="k",...) 
    if (length(x$n[which(!is.na(x$U.base))])<=3) points(x$n[which(!is.na(x$U.base))],x$U.base[which(!is.na(x$U.base))],col=1)
    else lines(x$n[which(!is.na(x$U.base))],x$U.base[which(!is.na(x$U.base))],col=1)

    arrows(x0=n.target-ESS, y0=Ub(n.target), x1 = n.target, y1 = Ub(n.target),col=2,length=.15)
    text(n.target-ESS/2,Ub(n.target),paste0("ECSS=",round(ESS,1)),cex=1,font=3,pos=3)
    
    abline(v=n.target,col=1,lty="dashed",lwd=2)
    abline(v=n.target-ESS,col=2)

    legend("topright",legend = c("base","info"),text.col = c(1,2))
    return(ESS)
}
