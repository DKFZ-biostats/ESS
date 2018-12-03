# dlink.mod <-
# function (object, value) {
#     if (RBesT:::is.dlink(value))
#         trans <- value
#     else {
#         trans <- match.fun(value)()
#     }
#   #  assert_that(is.dlink(trans))   #uncommenting this is the only change #   assertthat::assert_that(RBesT:::is.dlink(trans))
#    attr(object, "link") <- trans
#     object
# }


my.log_sum_exp=function(x) {
    if (length(x) == 1) return(x)
  
    xmax <- which.max(x)
    if (is.finite(x[xmax])) return(log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax])
    
    return(log(sum(exp(x))))
}
 
my.postmix.normMix <-
function (priormix, data, n, m, se, ...) { #allow vector m
    dataPrec <- 1/se^2
    priorPrec <- 1/priormix[3, , drop = FALSE]^2
    postPrec <- priorPrec + dataPrec
    sigmaPred <- sqrt(priormix[3, , drop = FALSE]^2 + se^2)
    logpriormean=log(priormix[1, , drop = FALSE])
    weight= sapply(1:ncol(priormix), function(cc) dnorm(m, priormix[2,cc
        , drop = T], sigmaPred[cc], log = TRUE))
    
    lwmat <- sapply(1:ncol(priormix), function(cc)  logpriormean[cc]+weight[,cc]) #switch(is.vector(weight)+1,weight[,cc],weight[cc])
    mix1=apply(lwmat,1,function(lw) exp(lw -my.log_sum_exp(lw)))
    #priormix[1, ] <- exp(lw -RBesT:::log_sum_exp(lw))
    mix2=t(sapply(1:ncol(priormix), function(cc) (priormix[2,cc , drop = T] * priorPrec[cc] +
        m * dataPrec)/postPrec[cc]))
 #       priormix[2, ] <- (priormix[2, , drop = FALSE] * priorPrec +
  #      m * dataPrec)/postPrec

    postsd=1/sqrt(postPrec)
   # priormix[3, ] <- 1/sqrt(postPrec)
   res= lapply(1:length(m), function(mm){
     # a=priormix
     # a[1,]=mix1[,mm]
     # a[2,]=mix2[,mm]
     # a[3,]=postsd
     # class(a) <- c("normMix", "mix")
     #  a
     rbind(mix1[,mm],mix2[,mm],postsd)
    })
    res
}

# my.postmix.normMix <-
# function (priormix, data, n, m, se, ...) {
#     if (!missing(data)) {
#         m <- mean(data)
#         n <- length(data)
#         se <- sd(data)/sqrt(n)
#     }
#     else {
#         # if (missing(m) & (missing(n) | missing(se)))
#         #     stop("Either raw data or summary data (m and se) must be given.")
#         if (!missing(se) & !missing(n)) {
#             sigma(priormix) <- se * sqrt(n)
#             # message(paste0("Updating reference scale to ", sigma(priormix),
#             #     ".\nIt is recommended to use the sigma command instead.\nSee ?sigma or ?mixnorm."))
#         }
#         if (missing(se) & !missing(n)) {
#      #       message("Using default prior reference scale ", sigma(priormix))
#             se <- sigma(priormix)/sqrt(n)
#         }
#     }
#     dataPrec <- 1/se^2
#     priorPrec <- 1/priormix[3, , drop = FALSE]^2
#     postPrec <- priorPrec + dataPrec
#     sigmaPred <- sqrt(priormix[3, , drop = FALSE]^2 + se^2)
#     lw <- log(priormix[1, , drop = FALSE]) + dnorm(m, priormix[2,
#         , drop = FALSE], sigmaPred, log = TRUE)
#     priormix[1, ] <- exp(lw -RBesT:::log_sum_exp(lw))
#     priormix[2, ] <- (priormix[2, , drop = FALSE] * priorPrec +
#         m * dataPrec)/postPrec
#     priormix[3, ] <- 1/sqrt(postPrec)
#     class(priormix) <- c("normMix", "mix")
#     priormix
# }
