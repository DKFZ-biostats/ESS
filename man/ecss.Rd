\name{ecss}
\alias{ecss}
\alias{ecss.mix}
\alias{ecss.mixture.prior}
\alias{ecss.ecss}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Effective Current Sample Size (ECSS)
}
\description{
Compute Effective Current Sample Size (ECSS)
}
\usage{
\method{ecss}{mix}(prior,
    data, n, r, m, sigma, se, true.mean,
    n.target, 
    min.ecss, 
    prior.base, D = MSE, 
    by, grid.length = 50, grid.prec = 1e-06, 
    cores = 1, integrate = FALSE, subdivisions = 100L, ...)

\method{ecss}{mixture.prior}(prior, ...)

\method{ecss}{ecss}(prior, n.target, tol = 1e-06, tol2 = 0.001, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{An RBesT betaMix or normMix object, a powerprior object created by as.powerprior, a StudyPrior mixture.prior object or an ecss object created by \code{ecss}
%%     ~~Describe \code{mix} here~~
}
  \item{data}{
individual data as in \code{\link[RBesT]{postmix}}. If the individual data is not given, then summary data has to be provided}
  \item{n}{
sample size}
  \item{r}{
in case of binary outcome, number of successes
}
  \item{m}{
in case of normal outcome, data mean
}
  \item{sigma}{
in case of normal outcome, sample standard deviation
}
  \item{se}{
in case of normal outcome, sample standard error
}
  \item{true.mean}{
  replaces m or r/n if mean is assumed to be known. In this case the truth is assumed to equal true.mean instead of the posterior mean given the data
}
  \item{n.target}{
  The sample size of interest at which the ECSS is evaluated. Can be a vector.
%%     ~~Describe \code{n.target} here~~
}
  \item{min.ecss}{
minimal ECSS of interest (negative)
}
  \item{prior.base}{
baseline prior
}
  \item{D}{
A function that measures informatives, e.g. \code{\link{MSE}}. % or \code{\link{BMSE}}.
}
  \item{by}{
  Defaults to 5 for normal outcome and to 1 for binary outcome. In case of binary outcome, by=1 recommended.
%%     ~~Describe \code{by} here~~
}
  \item{grid.length}{
in case of normal outcome, the grid length for integration
}
  \item{grid.prec}{
in case of normal outcome, small for a wide grid
}
  \item{cores}{
number of parallel cores used in mclapply
}
  \item{integrate}{
use \code{\link[stats]{integrate}}? Not recommended
}
  \item{subdivisions}{
in case \code{\link[stats]{integrate}} is used
}

%  \item{type}{
%don't change
%}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{tol2}{
%%     ~~Describe \code{tol2} here~~
}

  \item{\dots}{
  possible arguments passed to function provided in \code{D}
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Prints ECSS at n.target, returns an ecss object which can be updated using \link{ecss.ecss}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Wiesenfarth&Calderazzo (2019)
}
\author{
Manuel Wiesenfarth
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{ehss}
}
\examples{
######################
# Normal Outcome
  # data SD
    sigma=1
  # baseline
    rob=c(0,10)
    vague <-mixnorm(vague=c(1, rob), sigma=sigma)
  # prior with nominal EHSS=50
    inf=c(0,1/sqrt(50))
    info <-mixnorm(informative=c(1, inf), sigma=sigma)
  # robust mixture
    mix50 <-mixnorm(informative=c(.5, inf),vague=c(.5, rob), sigma=sigma)
  
  

  m=.2 #data mean
  n=100 # sample size
  
  ecss(info,m=m,sigma=sigma,n=n,
       n.target=100,min.ecss = -100,prior.base = vague)
  ecss(as.powerprior(info),m=m,sigma=sigma,n=n,
       n.target=100,min.ecss = -100,prior.base = vague)
  ecss(mix50,m=m,sigma=sigma,n=n,
       n.target=100,min.ecss = -100,prior.base = vague)
   
       
#############
# binary outcome
  unif  <- mixbeta(               c(1.0, 1, 1))
  beta  <- mixbeta(c(1.0, 4, 16)              )
  mix50 <- mixbeta(c(0.5, 4, 16), c(0.5, 1, 1))

  r=12
  n=20
  ecss(beta,r=r,n=n,
         n.target=20,min.ecss = -40,prior.base = unif,by=1)
  ecss(mix50,r=r,n=n,
         n.target=20,min.ecss = -10,prior.base = unif,by=1)
  ecss(as.powerprior(beta,p.prior.a=1,p.prior.b=1),r=r,n=n,
         n.target=20,min.ecss = -10,prior.base = unif,by=1)


###############
#library(StudyPrior)
#pars <- matrix(c(2.4,1,3.7,1),ncol=2)
#
##A mixture of 0.5 Be(2.4,3.7) + 0.5 Be(1,1) 
#mymix <- create.mixture.prior(type="beta", pars = pars, weights=c(0.5,0.5))
#ecss(mymix,r=r,n=n,
#         n.target=20,min.ecss = -10,prior.base = unif,by=1)



}
