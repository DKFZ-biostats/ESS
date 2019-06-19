\name{as.commensurateEB}
\alias{as.commensurateEB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert RBesT mix object into a commensurateEB object
}
\description{
Convert RBesT mix object into commensurateEB object allowing the computation of empirical Bayes commensurate priors.
}
\usage{
as.commensurateEB(prior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{An RBesT betaMix or normMix mixture object with a single mixture component
%%     ~~Describe \code{prior} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A commensurateEB object which can be used in \code{\link{ecss}}. %, e.g. in \code{\link{postmix.commensurateEB}} to compute posterior. 
Note that for the normal-normal model, the empirical Bayes power and commensurate prior coincide (see Wiesenfarth and Calderazzo, 2019).
%FOR BETA-BINOMIAL MODEL BETA APPROXIMATION IS USED

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hobbs, B. P., Sargent, D. J., and Carlin, B. P. (2012).  Commensurate priors for incorporating
historical information in clinical trials using general and generalized linear models.
Bayesian
analysis (Online)
7,
639.

Wiesenfarth, M., Calderazzo, S. (2019). Quantification of Prior Impact in Terms of Effective Current Sample Size. Submitted.



}
\author{
Manuel Wiesenfarth
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\link{as.powerprior}
}
\examples{
\dontrun{
#############
# binary outcome
  unif  <- mixbeta(               c(1.0, 1, 1))
  beta  <- mixbeta(c(1.0, 4, 16)              )

  r=7
  n=20
  ecss(as.commensurateEB(beta),r=r,n=n,
         n.target=40,min.ecss = 0,prior.base = unif,by=1)
}
}


