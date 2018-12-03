\name{adaptiveDesign_binomial}
\alias{adaptiveDesign_binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate adaptive design with adaptive adjustment of the randomiziation ratio according to ESS for binomial outcome
}
\description{
Simulate adaptive design with adaptive adjustment of the randomiziation ratio according to ESS for binomial outcome as in Schmidli at al (2014).
}
\usage{
adaptiveDesign_binomial(ctl.prior, treat.prior, 
    N1, Ntarget, Nmin, M, 
    pc, pt, 
    discard.prior = TRUE, 
    vague = mixbeta(c(1, 1, 1)), 
    ess = "ecss", 
    ehss.method = "mix.moment", substractEHSSofVague=TRUE,
    min.ecss, D=MSE,
    decision)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctl.prior}{
RBesT betaMix object (or powerprior object created by \code{\link{as.powerprior}}) as prior for the control group}
  \item{treat.prior}{
RBesT betaMix object (or powerprior object created by \code{\link{as.powerprior}}) as prior for the treatment group
}

  \item{N1}{
Sample size in each group at interim
}
  \item{Ntarget}{
Target sample size in control group
}
  \item{Nmin}{
Minimal number of samples in control group after interim analysis
}
  \item{M}{
Final sample size in treatment group
}
  \item{pc}{
True control rate
}
  \item{pt}{
True treatment rate
}
  \item{discard.prior}{
Replace prior by vague prior if ESS<0?
}
  \item{vague}{
RBesT betaMix object serving as baseline vague prior
}
  \item{ess}{
either "ecss" or "ehss" for effective current or historical sample size, respectively.
}
  \item{ehss.method}{
if ess=="ehss". Specify version of EHSS as in \code{\link{ehss}}.
}
  \item{substractEHSSofVague}{
}

  \item{min.ecss}{
if ess=="ecss". Minimal ECSS of interest (negative)
}
  \item{D}{
  A function that measures informatives, e.g. \code{\link{MSE}}
}
  \item{decision}{
function created by \code{\link[RBesT]{decision2S}}.
}
}


%\details{
%}
%\value{
%}


\references{
Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., and Neuenschwan-
der, B. (2014). Robust meta-analytic-predictive priors in clinical trials with historical control
information. Biometrics, 70(4):1023-103
}
\author{
Manuel Wiesenfarth
}
%\note{
%}


\seealso{
\code{vignette("robustMAP",package="RBesT")}
}

%\examples{
%}

