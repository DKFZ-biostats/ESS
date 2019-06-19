\name{adaptiveDesign_normal}
\alias{adaptiveDesign_normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate adaptive design where control sample size is adjusted according to ESS for normal outcome
}
\description{
Simulate adaptive design where control sample size is adapted according to prior effective sample size for normal outcome similar to Schmidli at al (2014). 
}
\usage{
adaptiveDesign_normal(ctl.prior, treat.prior, 
    N1, Ntarget, Nmin, M, 
    muc, mut, sc, st, sc.known=TRUE, st.known=TRUE,
    discard.prior = TRUE, 
    vague = mixnorm(c(1, 0, 10)), 
    ess = "ecss", ehss.method = "mix.moment", min.ecss, D=MSE,
    decision, 
    nsim = 100, cores = 1, seed = 123, progress="text")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctl.prior}{
RBesT normMix object (or powerprior object created by \code{\link{as.powerprior}}) as prior for the control group}
  \item{treat.prior}{
RBesT normMix object (or powerprior object created by \code{\link{as.powerprior}}) as prior for the treatment group
}
  \item{N1}{
Sample size in each group at interim
}
  \item{Ntarget}{
Target sample size in control group
}
  \item{Nmin}{
Minimum number of samples in control group after interim analysis
}
  \item{M}{
Final sample size in treatment group
}

  \item{muc}{
True control mean}
  \item{mut}{
True mean in treatment group
}
  \item{sc}{
standard deviation in control group
}
  \item{st}{
standard deviation in treatment group
}
  \item{sc.known}{
  logical. If \code{TRUE}, assume \code{sc} to be known, otherwise replace by empirical standard deviation.
}
  \item{st.known}{
  logical. If \code{TRUE}, assume \code{st} to be known, otherwise replace by empirical standard deviation.
}
  \item{discard.prior}{
Replace prior by vague prior if ESS<0?
}
  \item{vague}{
RBesT normMix object (single component mixture prior) serving as baseline vague prior
}
  \item{ess}{
either "ecss" or "ehss" for effective current or historical sample size, respectively.
}
  \item{ehss.method}{
if ess=="ehss". Specify version of EHSS as in \code{\link{ehss}}.
}
  \item{min.ecss}{
if ess=="ecss". Minimal ECSS of interest (negative). A large absolute value of min.ecss is computational expensive, could be set to -1 if \code{discard.prior=TRUE} and no interest in the ECSS estimate itself.
}
  \item{D}{
  A function that measures informatives, e.g. \code{\link{MSE}} or user-specified function
}
  \item{decision}{
function created by \code{\link[RBesT]{decision2S}}.
}
  \item{nsim}{
Number of Monte Carlo iterations
}
  \item{cores}{
number of parallel cores used in mclapply
}
  \item{seed}{
random seed
}
  \item{progress}{
name of the progress bar to use (see \code{\link[plyr]{create_progress_bar}})
}

}
\details{
The traditional approach to prior effective sample size (prior ESS) is aimed at quantifying prior informativeness, but is not aimed at detecting potential prior-data conflict.

The ECSS computes the prior effective sample size in terms of samples from the current data model (i.e., samples with characteristics consistent with the current trial). Under extreme prior-data conflict, the prior may account for a negative number of samples, showing that information is subtracted, rather than added, by the elicited prior. The ECSS quantifies the number of current samples to be added or subtracted to the likelihood in order to obtain a posterior inference equivalent to that of a baseline prior model (e.g. in terms of mean squared error, MSE). Fur further details, see Wiesenfarth and Calderazzo (2019).

Standard approach uses effective historical sample size (\code{ess="ehss"}), while Wiesenfarth and Calderazzo (2019) use the effective current sample size (\code{ess="ecss"}).
When the ECSS is negative, the design provides the option of discarding the prior (discard.prior=TRUE). 

Extensive documentation is given in the vignette.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., and Neuenschwan-
der, B. (2014). Robust meta-analytic-predictive priors in clinical trials with historical control
information. Biometrics, 70(4):1023-103

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
\code{vignette("robustMAP",package="RBesT")}
}

\examples{
# see
# vignette("vignetteDesign", package = "ESS")
}
