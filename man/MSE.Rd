\name{MSE}
\alias{MSE}
\alias{MSEmix.betaMix}
\alias{MSEmix.normMix}
\alias{MSEuni.betaMix}
\alias{MSEuni.normMix}
\alias{MSEpp.betaMix}
\alias{MSEpp.normMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute frequentist MSE
}
\description{
Compute frequentist MSE for use as argument in \code{\link{ecss}}
}
\usage{
MSE(prior, ...)
\method{MSEmix}{betaMix}(prior, y, k, true.mean, ...)
\method{MSEmix}{normMix}(prior, y, k, true.mean, ...)
\method{MSEuni}{betaMix}(prior, y, k, true.mean, ...)
\method{MSEuni}{normMix}(prior, y, k, true.mean, ...)
\method{MSEpp}{betaMix}(prior, y, k, true.mean, ...)
\method{MSEpp}{normMix}(prior, y, k, true.mean, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{An RBesT betaMix or normMix object, or a powerprior object created by \link{as.powerprior}.
%%     ~~Describe \code{PRIOR} here~~
}
  \item{y}{
vector of quantiles
}
  \item{k}{
number of observations
}
  \item{true.mean}{
true mean or rate
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
%\examples{

%}

