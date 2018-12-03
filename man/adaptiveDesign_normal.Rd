\name{adaptiveDesign_normal}
\alias{adaptiveDesign_normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate adaptive design with adaptive adjustment of the randomiziation ratio according to ESS for normal outcome
}
\description{
Simulate adaptive design with adaptive adjustment of the randomiziation ratio according to ESS for normal outcome similar to Schmidli at al (2014).
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
Minimal number of samples in control group after interim analysis
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
RBesT normMix object serving as baseline vague prior
}
  \item{ess}{
either "ecss" or "ehss" for effective current or historical sample size, respectively.
}
  \item{ehss.method}{
if ess=="ehss". Specify version of EHSS as in \code{\link{ehss}}.
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
name of the progress bar to use
}

}
%\details{
%}
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
}
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

