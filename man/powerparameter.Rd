\name{powerparameter}
\alias{powerparameter}
\alias{powerparameter.betaMix}
\alias{powerparameter.normMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes empirical Bayes esitmate of power parameter
%%  ~~function to do ... ~~
}
\description{
Computes empirical Bayes esitmate of power parameter (Gravestock&Held, 2017).
}
\usage{
\method{powerparameter}{betaMix}(prior, n, r, p.prior.a = 1, p.prior.b = 1, ...)
\method{powerparameter}{normMix}(prior, n, m, sigma, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{An RBesT mixture object with a single mixture component
%%     ~~Describe \code{prior} here~~
}
  \item{n}{
sample size}
  \item{r}{
Number of successes
}
  \item{m}{
Sample mean
}
  \item{sigma}{
Sample standard deviation
}
  \item{p.prior.a}{
in case of binary outcome, shape1 parameter of initial beta prior for successes. 
}
  \item{p.prior.b}{
in case of binary outcome, shape2 parameter of initial beta prior for successes. 
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

\seealso{ package StudyPrior, 
\link{postmix.powerprior}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Normal Outcome
  # data SD
    sigma=1
  # prior with nominal EHSS=50
    info <-mixnorm(informative=c(1, 0,1/sqrt(50)), sigma=sigma)
  
  n=10 # sample size

  m=0 #data mean equal to prior mean
  powerparameter(info,m=m,n=n,sigma=sigma)
  m=2 #prior-data conflict
  powerparameter(info,m=m,n=n,sigma=sigma)

}

