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
Compute Effective Current Sample Size (ECSS) which - roughly - quantifies the number of current samples to be added or subtracted to the likelihood in order to obtain a posterior inference equivalent to that of a baseline prior model (e.g. in terms of mean squared error, MSE). 
}
\usage{
\method{ecss}{mix}(prior,
    data, n, r, m, sigma, se, true.mean,
    n.target, 
    min.ecss, 
    prior.base, D = MSE, 
    by, grid.length = 50, min.q, 
    cores = 1, integrate = FALSE, subdivisions = 100L, ...)

\method{ecss}{mixture.prior}(prior, ...)

\method{ecss}{ecss}(prior, n.target, tol = 1e-06, tol2 = 0.001, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{An \code{RBesT betaMix} or \code{normMix} object, a powerprior object created by \code{\link{as.powerprior}}, a \code{StudyPrior} mixture.prior object or an \code{ecss} object created by \code{ecss}
%%     ~~Describe \code{mix} here~~
}
  \item{data}{
if applied to given data, individual data as in \code{\link[RBesT]{postmix}}. If the individual data is not given, then summary data has to be provided}
  \item{n}{
if applied to given data, sample size}
  \item{r}{
in case of binary outcome, if applied to given data, number of successes
}
  \item{m}{
in case of normal outcome, if applied to given data, data mean
}
  \item{sigma}{
in case of normal outcome, sample standard deviation
}
  \item{se}{
in case of normal outcome, sample standard error
}
  \item{true.mean}{the true mean in case of prospective use at planning stage.
 % replaces m or r/n if mean is assumed to be known. 
  If provided, the truth is assumed to equal true.mean instead of the posterior mean given the data
}
  \item{n.target}{
  The sample size of interest at which the ECSS is evaluated. Can be a vector.
%%     ~~Describe \code{n.target} here~~
}
  \item{min.ecss}{
minimal ECSS of interest (usually negative)
}
  \item{prior.base}{RBesT betaMix or normMix object (single mixture component) serving as baseline prior. Usually flat prior and uniform prior in case of normal and beta prior of interest, respectively.
}
  \item{D}{
A function that measures informativeness, e.g. \code{\link{MSE}} or user-specified function. % or \code{\link{BMSE}}.
}
  \item{by}{stepsize of effective sample sizes where \code{D} is evaluated, interpolation in between.
  Defaults to 5 for normal outcome and to 1 for binary outcome. In case of binary outcome, by=1 recommended.
%%     ~~Describe \code{by} here~~
}
  \item{grid.length}{
number of elements on grid for integration, with bounds determined by \code{min.q}
}
  \item{min.q}{
lower quantile of normal or binomial likelihood to specify width of grid for integration (1-min.q upper quantile), small for a wide grid. Defaults to 1e-6 for \code{normMix} objects and 0 for \code{betaMix} objects.
}
  \item{cores}{
number of parallel cores used in mclapply
}
  \item{integrate}{
use \code{\link[stats]{integrate}}? Not recommended
}
  \item{subdivisions}{
passed to \code{\link[stats]{integrate}} if \code{integrate=FALSE}
}
  \item{tol}{
  not to be changed by user
}
  \item{tol2}{
  only change if result is NA, and use \code{\link{plot.ecss}} for diagnostics then
}
  \item{\dots}{
  possible arguments passed to function provided in \code{D}
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Let \eqn{\pi} be the prior of interest specified by \code{prior} with mean \eqn{\theta_\pi}, \eqn{\pi_b} a baseline prior (an objective or reference prior) specified by \code{prior.base} and \eqn{f_n(y_{1:n} | \theta_0)} be the data distribution. 

The ECSS at target sample size \eqn{k} (specified by \code{n.target}) is defined as the sample size \eqn{m} which minimizes
\deqn{ECSS=argmin_{m} | D^{\theta_0}_{MSE}({\pi}(\theta|y_{1:(k-m)})) -  D^{\theta_0}_{MSE}(\pi_b(\theta|y_{1:k})) |}

where \eqn{D^{\theta_0}_{MSE}} is the mean squared error measure induced by the posterior mean estimate \eqn{\mbox{E}_\pi(\theta|y_{1:k})},
\deqn{D^{\theta_0}_{MSE}({\pi}(\theta|y_{1:k}))=E_{y|\theta_0 } [ E_\pi(\theta|y_{1:k})-\theta_0 ]^2}
or a different target measure specified by \code{D}.


The true parameter value \eqn{\theta_0} is either assumed to be known and specified by \code{true.mean} (useful for prospective quantification of the prior impact at planning stage) or replaced by the posterior mean under the baseline prior specification given specified data (useful for sensitivity analyses or adaptive designs).

See Wiesenfarth and Calderazzo (2019) for details.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
Prints ECSS at n.target, returns an ecss object which can be updated using \link{ecss.ecss}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Wiesenfarth, M., Calderazzo, S. (2019). Quantification of Prior Impact in Terms of Effective Current Sample Size. Submitted.

}
\author{
Manuel Wiesenfarth
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{ess}, \link{ehss}
}
\examples{
######################
# Normal Outcome
  # standard deviation
    sigma=1
  # baseline
    rob=c(0,10)
    vague <-mixnorm(vague=c(1, rob), sigma=sigma)
  # prior with nominal prior ESS=50
    inf=c(0,1/sqrt(50))
    info <-mixnorm(informative=c(1, inf), sigma=sigma)
  # robust mixture
    mix50 <-mixnorm(informative=c(.5, inf),vague=c(.5, rob), sigma=sigma)
  
  
  # before seeing data
    true.mean=.2 # hypothetic true value
    ecss(info,true.mean=true.mean,sigma=sigma,
         n.target=100,min.ecss = -100,prior.base = vague)
    ecss(as.powerprior(info),true.mean=true.mean,sigma=sigma,
         n.target=100,min.ecss = -100,prior.base = vague)
    ecss(mix50,true.mean=true.mean,sigma=sigma,
         n.target=100,min.ecss = -100,prior.base = vague)
  
  # for observed data
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

  # before seeing data
    true.mean=10/20
    ecss(beta,true.mean=true.mean,
           n.target=20,min.ecss = -40,prior.base = unif,by=1)
    ecss(mix50,true.mean=true.mean,
           n.target=20,min.ecss = -10,prior.base = unif,by=1)
    ecss(as.powerprior(beta,p.prior.a=1,p.prior.b=1),true.mean=true.mean,
           n.target=20,min.ecss = -10,prior.base = unif,by=1)


  # for observed data
    r=10
    n=20
    ecss(beta,r=r,n=n,
           n.target=20,min.ecss = -40,prior.base = unif,by=1)
    ecss(mix50,r=r,n=n,
           n.target=20,min.ecss = -10,prior.base = unif,by=1)
    ecss(as.powerprior(beta,p.prior.a=1,p.prior.b=1),r=r,n=n,
           n.target=20,min.ecss = -10,prior.base = unif,by=1)


\dontrun{
library(StudyPrior)
pars <- matrix(c(2.4,1,3.7,1),ncol=2)

#A mixture of 0.5 Be(2.4,3.7) + 0.5 Be(1,1) 
mymix <- create.mixture.prior(type="beta", pars = pars, weights=c(0.5,0.5))
ecss(mymix,r=r,n=n,
         n.target=20,min.ecss = -10,prior.base = unif,by=1)
}


}
