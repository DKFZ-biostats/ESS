\name{postmix.powerprior}
\alias{postmix.powerprior}
\alias{powerprior.betaMix}
\alias{powerprior.normMix}
\alias{powerprior}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Analysis for empirical Bayes powerprior
}

\description{
Calculates the posterior distribution given data \code{data} and prior \code{priormix}, where the prior is an empirical Bayes powerprior (Gravestock&Held, 2017). Supports beta-binomial and normal-normal model.
}
\usage{
\method{postmix}{powerprior}(priormix, data, n, r, m, se, p.prior.a, p.prior.b, ...)
\method{powerprior}{betaMix}(prior, n, r, p.prior.a, p.prior.b, ...)
\method{powerprior}{normMix}(prior, n, m, sigma, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{priormix}{
prior object created by \code{\link{as.powerprior}}
}
 \item{prior}{An RBesT mixture object with a single mixture component}
  \item{data}{
individual data as in \code{\link[RBesT]{postmix}}. If the individual data is not given, then summary data has to be provided
}
  \item{n}{
sample size}
  \item{r}{
number of successes
}
  \item{m}{
sample mean
}
  \item{se}{
sample standard error
}
\item{sigma}{standard deviation}
  \item{p.prior.a}{
  in case of binary outcome, shape1 parameter of initial beta prior for successes
}
  \item{p.prior.b}{
  in case of binary outcome, shape2 parameter of initial beta prior for successes
}
  \item{\dots}{
  currently not supported
}
}
\details{
Extends \code{\link[RBesT]{postmix}}. See its documentation for further details on arguments.
\code{powerprior} and \code{postmix.powerprior} are equivalent except that the former accepts a single component mixture object while the latter accepts a \code{powerprior} object created by \code{\link{as.powerprior}}. \code{powerprior} may be deprecated in the future.
%%  ~~ If necessary, more details than the description above ~~
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Gravestock, I. and Held, L. (2017). Adaptive power priors with empirical bayes for clinical trials.
Pharmaceutical statistics, 16(5):349-360.
}
\author{
Manuel Wiesenfarth
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{package StudyPrior
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
######################
# Normal Outcome
  # standard deviation
    sigma=1
  # prior with nominal prior ESS=50
    inf=c(0,1/sqrt(50))
    info <-mixnorm(informative=c(1, inf), sigma=sigma)

  n=10
  posterior=postmix(as.powerprior(info),n = n,m=1,se=sigma/sqrt(n))
  plot(posterior)
  pmix(posterior, 1.5, lower.tail = FALSE)
}
