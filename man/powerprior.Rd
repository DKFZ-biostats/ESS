\name{powerprior}
\alias{powerprior.betaMix}
\alias{powerprior.normMix}
\alias{powerprior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Analysis for empirical Bayes powerprior
}
\description{
Calculates the posterior distribution for data data given a prior priormix, where the prior is a empirical Bayes power prior (Gravestock&Held, 2017). Does the same as \link{postmix.powerprior}. Maybe deprecated in the future.
}
\usage{
\method{powerprior}{betaMix}(prior, n, r, p.prior.a, p.prior.b, ...)
\method{powerprior}{normMix}(prior, n, m, sigma, ...)
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

\seealso{package StudyPrior, 
\link{postmix.powerprior}
}
%\examples{

%}
