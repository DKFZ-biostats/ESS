\name{ehss}
\alias{ehss}
\alias{ehss.normMix}
\alias{ehss.betaMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Effective Historical Sample Size (EHSS)
}
\description{
Compute Effective Historical Sample Size (EHSS)
}
\usage{
\method{ehss}{normMix}(prior, method = c("mix.moment", "moment", "morita"), ...)
\method{ehss}{betaMix}(prior, method = c("mix.moment", "moment", "morita"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{
An RBesT mixture object
}
  \item{method}{
  Selects the used method. Can be either mix.moment, moment or morita.
%%     ~~Describe \code{method} here~~
}
   \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Methods "moment" and "morita" calculates the moments based approach and the method by Morita et al. (2008) using implementation in \code{\link[RBesT]{ess}}.
"mix.moment" computes the EHSS for each mixture component and weights them according to the mixture weights as e.g. in Gravestock&Held (2018).
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
Gravestock, I. and Held, L. (2018).  Power priors based on multiple historical studies for binary
outcomes. Biometrical Journal, 0(0)

Morita S, Thall PF, Mueller P. Determining the effective sample size of a parametric prior. Biometrics 2008;64(2):595-602.
}
\author{
Manuel Wiesenfarth
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[RBesT]{ess}, \link{ecss}
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
  
  
  ehss(info,method="moment")
  
  ehss(mix50,method="morita")
  ehss(mix50,method="moment")
  ehss(mix50,method="mix.moment")
  
  
  m=.2 #data mean
  n=100 # sample size
  
  ehss(powerprior(info,m=m,n=n,sigma=sigma)) 
  ehss(postmix(as.powerprior(info),m=m,n=n,se=sigma/sqrt(n)))-n

}
