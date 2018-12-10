\name{plotECSS}
\alias{plotECSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot ECSS 
}
\description{
Plot ECSS for a list of priors on a grid
}
\usage{
plotECSS(priorlist, grid, n.target, min.ecss, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{priorlist}{
A named list of priors
}
  \item{grid}{
%%     ~~Describe \code{grid} here~~
}
  \item{n.target}{
%%     ~~Describe \code{n.target} here~~
}
  \item{min.ecss}{
%%     ~~Describe \code{min.ecss} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{plot}{A "ggplot" object}
  \item{data }{A data frame with computed values}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Manuel Wiesenfarth
}

\examples{
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

plotECSS(priorlist=list(mix=mix50,info=info), grid=((-8:8)/10), n.target=200, min.ecss=-200,sigma=sigma,prior.base=vague)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
