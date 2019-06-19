\name{as.powerprior}
\alias{as.powerprior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert RBesT mix object into powerprior object
}
\description{
Convert RBesT mix object into powerprior object allowing the computation of empirical Bayes power priors (Gravestock&Held, 2017).
}
\usage{
as.powerprior(prior, p.prior.a = 1, p.prior.b = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{An RBesT betaMix or normMix mixture object with a single mixture component
%%     ~~Describe \code{prior} here~~
}
  \item{p.prior.a}{
  in case of binary outcome, shape1 parameter of initial beta prior for successes
}
  \item{p.prior.b}{
  in case of binary outcome, shape2 parameter of initial beta prior for successes
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A powerprior object which can be used, e.g. in \code{\link{postmix.powerprior}} to compute posterior.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
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
    info <-mixnorm(informative=c(1, 0,1/sqrt(50)), sigma=sigma)
  
  n=10 # sample size

  m=0 #data mean equal to prior mean
  postmix(info,m=m,se=sigma/sqrt(n))
  postmix(as.powerprior(info),m=m,se=sigma/sqrt(n))
  m=2 #prior-data conflict
  postmix(info,m=m,se=sigma/sqrt(n))
  postmix(as.powerprior(info),m=m,se=sigma/sqrt(n)) #variance inflation

}


