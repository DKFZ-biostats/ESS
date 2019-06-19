\name{plotECSS}
\alias{plotECSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot effective current sample size (ECSS) on a grid 
}
\description{
Plot effective current sample size (ECSS) on a grid for a list of priors.
}
\usage{
plotECSS(priorlist, grid, n.target, min.ecss, progress = "text", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{priorlist}{
A named list of priors. Elements can be \code{RBesT betaMix} or \code{normMix} objects, powerprior objects created by \code{\link{as.powerprior}} or \code{StudyPrior} mixture.prior objects.
}
  \item{grid}{
  The grid of true data means (see details in \code{link{ecss.mix}}).
%%     ~~Describe \code{grid} here~~
}
  \item{n.target}{
  The sample size of interest at which the ECSS is evaluated
%%     ~~Describe \code{n.target} here~~
}
  \item{min.ecss}{
minimal ECSS of interest (usually negative)
}
  \item{progress}{
name of the progress bar to use, see \code{\link[plyr]{create_progress_bar}} 
}
  \item{\dots}{
  Further arguments to be passed to \code{\link{ecss}}.
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Plot effective current sample size (ECSS) on a grid for a list of priors. This is useful to visualize the impact of a prior for a range of true values at planning stage, i.e. in prospective use before seeing the data. See also vignette and Wiesenfarth and Calderazzo (2019).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{plot}{A "ggplot" object}
  \item{data }{A data frame with computed values}
%% ...
}
\references{
Wiesenfarth, M., Calderazzo, S. (2019). Quantification of Prior Impact in Terms of Effective Current Sample Size. Submitted.

}
\author{
Manuel Wiesenfarth
}

\examples{
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

plotECSS(priorlist=list(mix=mix50,info=info), grid=((-8:8)/10), n.target=200, 
          min.ecss=-200,sigma=sigma,prior.base=vague)

# See also vignette for more extensive example.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
