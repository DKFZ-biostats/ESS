\name{plot.ecss}
\alias{plot.ecss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot ecss object
}
\description{
Plot an ecss object created by \code{ecss} for diagnostics. 
}
\usage{
\method{plot}{ecss}(x, n.target, 
    tol = 1e-06, tol2 = 0.001,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An \code{ecss} object
%%     ~~Describe \code{x} here~~
}
  \item{n.target}{The sample size of interest at which the ECSS is evaluated.
%%     ~~Describe \code{k} here~~
}
%  \item{type}{
%}
  \item{tol}{not to be changed by user
%%     ~~Describe \code{tol} here~~
}
  \item{tol2}{only change if result is NA
%%     ~~Describe \code{tol2} here~~
}
  \item{\dots}{arguments passed to \code{plot.default}
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Takes an \code{ecss} object and draws \code{D} (e.g. MSE if specified when calling  \code{\link{ecss}}) for the prior of interest and baseline prior which can be useful to understand whether small changes could have changed estimation of the effective current sample size.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Manuel Wiesenfarth
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  # standard deviation
    sigma=1
  # baseline
    rob=c(0,10)
    vague <-mixnorm(vague=c(1, rob), sigma=sigma)
  # prior with nominal prior ESS=50
    inf=c(0,1/sqrt(50))
    info <-mixnorm(informative=c(1, inf), sigma=sigma)

  m=.2 #data mean
  n=100 # sample size
  
  e=ecss(info,m=m,sigma=sigma,n=n,
       n.target=1:100,min.ecss = -100,prior.base = vague)
  plot(e,n.target=100)

}
