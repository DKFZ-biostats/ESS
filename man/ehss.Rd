\name{ehss}
\alias{ehss}
\alias{ehss.normMix}
\alias{ehss.betaMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Effective Historical Sample Size (EHSS)
}
\description{
Compute Effective Historical Sample Size (EHSS). This is the prior effective sample size applied to the posterior subtracting the data sample size, see also Wiesenfarth and Calderazzo (2019). 
}
\usage{
\method{ehss}{normMix}(prior, data, n, m, se, method = c("mix.moment", "moment", "morita"), ...)
\method{ehss}{betaMix}(prior, data, n, r, method = c("mix.moment", "moment", "morita"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{
An RBesT mixture object
}
  \item{data}{
individual data as in \code{\link[RBesT]{postmix}}. If the individual data is not given, then summary data has to be provided
}
  \item{n}{
sample size}
  \item{r}{
number of successes
}
  \item{m}{
sample mean
}
  \item{se}{
sample standard error
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
Simply applies \code{\link{ess}} to the posterior and subtracts the sample size.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Wiesenfarth, M., Calderazzo, S. (2019). Quantification of Prior Impact in Terms of Effective Current Sample Size. Submitted.

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

  # standard deviation
    sigma=1
  # baseline
    rob=c(0,10)
    vague <-mixnorm(vague=c(1, rob), sigma=sigma)
  # prior with nominal EHSS=50
    inf=c(0,1/sqrt(50))
    info <-mixnorm(informative=c(1, inf), sigma=sigma)
  # robust mixture
    mix50 <-mixnorm(informative=c(.5, inf),vague=c(.5, rob), sigma=sigma)
  
  

  m=.2 #data mean
  n=100 # sample size
  
  ehss(as.powerprior(info),m=m,n=n,se=sigma/sqrt(n))
  
  ehss(mix50,m=m,n=n,se=sigma/sqrt(n),method="morita")
  ehss(mix50,m=m,n=n,se=sigma/sqrt(n),method="moment")
  ehss(mix50,m=m,n=n,se=sigma/sqrt(n),method="mix.moment")

}
