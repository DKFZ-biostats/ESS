\name{ess}
\alias{ess}
\alias{ess.normMix}
\alias{ess.betaMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Prior Effective Sample Size (prior ESS)
}
\description{
Compute prior effective sample size (prior ESS). This
this comprises different classical approches to ESS quantification, as e.g. the  approach of Morita et al (2008), see also Wiesenfarth and Calderazzo (2019).
}
\usage{
\method{ess}{normMix}(prior, method = c("mix.moment", "moment", "morita"), ...)
\method{ess}{betaMix}(prior, method = c("mix.moment", "moment", "morita"), ...)

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
The prior ESS can quantify the strength of prior information, but it is not designed to detect and/or quantify prior-data conflict. In a normal prior-normal likelihood model, for example, all priors sharing the same variance have the same prior ESS, but their impact in terms of posterior bias (and thus mean squared error, MSE) can vary substantially depending on their location with respect to the observed data. There is no consensus to date on how to compute a prior ESS for mixture priors and different options in \code{method} can provide considerably different results for the same prior, see also Wiesenfarth and Calderazzo (2019).

The prior is viewed as arising from a historical trial in which a certain number of observations were collected. Bayes' rule allows to incorporate such observations sequentially in the analysis, by considering the posterior distribution obtained in the historical trial as the prior distribution for the current trial. Well-known analytical results are available for conjugate priors of exponential family likelihoods, where such measure can be easily derived as a function of the prior parameters. Morita et al (2008) proposed an algorithmic extension applicable to arbitrary prior-likelihood settings, however, no consensus exists on its computation in case of mixture priors.   


Methods "moment" and "morita" calculates the moments based approach and the method by Morita et al. (2008), respectively, using implementation in \code{\link[RBesT]{ess}}, see documentation there.
"mix.moment" computes the prior ESS for each mixture component and weights them according to the mixture weights as e.g. in Gravestock&Held (2018). If the prior is a non-mixture (single component) beta or normal distribution, all these approaches yield the same results.
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

Morita S, Thall PF, Mueller P. (2008). Determining the effective sample size of a parametric prior. Biometrics, 64(2):595-602.

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
  # prior with nominal prior ESS=50
    inf=c(0,1/sqrt(50))
    info <-mixnorm(informative=c(1, inf), sigma=sigma)
  # robust mixture
    mix50 <-mixnorm(informative=c(.5, inf),vague=c(.5, rob), sigma=sigma)
  
  
  ess(info,method="moment")
  
  ess(mix50,method="morita")
  ess(mix50,method="moment")
  ess(mix50,method="mix.moment")
  
  

}
