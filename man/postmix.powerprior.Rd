\name{postmix.powerprior}
\alias{postmix.powerprior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Analysis for empirical Bayes powerprior
}

\description{
Calculates the posterior distribution for data data given a prior priormix, where the prior is a powerprior (Gravestock&Held, 2017). Extends \code{\link[RBesT]{postmix}}.
}
\usage{
\method{postmix}{powerprior}(priormix, data, n, r, m, se, p.prior.a, p.prior.b, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{priormix}{
prior object created by \code{\link{as.powerprior}}
}
  \item{data}{
individual data as in \code{\link[RBesT]{postmix}}. If the individual data is not given, then summary data has to be provided
}
  \item{n}{
sample size}
  \item{r}{
Number of successes
}
  \item{m}{
Sample mean
}
  \item{se}{
Sample standard error
}
  \item{p.prior.a}{
  in case of binary outcome, shape1 parameter of initial beta prior for successes
}
  \item{p.prior.b}{
  in case of binary outcome, shape2 parameter of initial beta prior for successes
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
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
  # data SD
    sigma=1
  # prior with nominal EHSS=50
    inf=c(0,1/sqrt(50))
    info <-mixnorm(informative=c(1, inf), sigma=sigma)

  n=10
  posterior=postmix(as.powerprior(info),n = n,m=1,se=sigma/sqrt(n))
  plot(posterior)
  pmix(posterior, 1.5, lower.tail = FALSE)
}
