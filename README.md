Quantification of prior impact in terms of effective current sample size
================

- [Visualization of the impact of priors in terms of effective sample
  size](#visualization-of-the-impact-of-priors-in-terms-of-effective-sample-size)
- [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

Computes prior effective and effective current sample size. The concept
of prior effective sample sizes (prior ESS) is convenient to quantify
and communicate informativeness of prior distributions as it equates the
information provided by a prior to a sample size. Prior information can
arise from historical observations, thus the traditional approach
identifies the ESS with such historical sample size (prior ESS).
However, this measure is independent from newly observed data, and thus
would not capture an actual \`\`loss of information’’ induced by the
prior in case of prior-data conflict. The effective current sample size
(ECSS) of a prior relates prior information to a number of (virtual)
samples from the current data model and describes the impact of the
prior capturing prior-data conflict. Supports mixture and empirical
Bayes power and commensurate priors. See Wiesenfarth and Calderazzo
(2019).

## Installation

<!-- To get the current released version from CRAN: -->
<!-- ```{r} -->
<!-- ## install ESS from CRAN -->
<!-- install.packages("ESS") -->
<!-- ## load ESS package -->
<!-- library(ESS) -->
<!-- ``` -->

To get the current development version from Github:

``` r
devtools::install_github("DKFZ-biostats/ESS")
```

``` r
## load ESS package
library(ESS)
#> Loading required package: RBesT
#> This is RBesT version 1.7.2 (released 2023-08-21, git-sha 0e02ce3)
#> 
#> ESS 0.5.0 loaded.
#> 
#> Attaching package: 'ESS'
#> The following object is masked from 'package:RBesT':
#> 
#>     ess
```

# Visualization of the impact of priors in terms of effective sample size

This section (same as `vignette("vignettePlanning", package = "ESS")`)
visualizes how `plotECSS()` can be used to quantify the impact of a
prior (including mixture priors and empirical Bayes power
priors/commensurate priors) in terms of effective current sample sizes
on a grid of true values of the data generating process.

The effective current sample size (ECSS) is defined as follows:

Let $\pi$ be the prior of interest (specified by `priorlist`) with mean
$\theta_\pi$, $\pi_b$ a baseline prior (an objective or reference prior)
(specified by `prior.base`) and $f_n(y_{1:n} | \theta_0)$ be the data
distribution.

The ECSS at target sample size $k$ (specified by `n.target`) is defined
as the sample size $m$ which minimizes
$$ECSS=argmin_{m} | D^{\theta_0}_{MSE}({\pi}(\theta|y_{1:(k-m)})) -  D^{\theta_0}_{MSE}(\pi_b(\theta|y_{1:k})) |$$

where $D^{\theta_0}_{MSE}$ is the mean squared error measure induced by
the posterior mean estimate $\mbox{E}_\pi(\theta|y_{1:k})$,
$$D^{\theta_0}_{MSE}({\pi}(\theta|y_{1:k}))=E_{y|\theta_0 } (E_\pi(\theta|y_{1:k})-\theta_0 )^2$$
or a different target measure specified by `D`. The true parameter value
$\theta_0$ (specified by `grid`) is assumed to be known.

## Normal Outcome

``` r
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
  # emprirical Bayes power prior / emprirical Bayes commensurate prior
    pp <-as.powerprior(info)

    plotECSS(priorlist=list(informative=info,mixture=mix50,powerprior=pp), 
             grid=((0:100)/100), n.target=200, min.ecss=-150, sigma=sigma, 
             prior.base=vague, progress="none", D=MSE)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Binary outcome

``` r
# uniform baseline prior
  rob=c(1,1)
  vague <- mixbeta(rob=c(1.0, rob))
# prior with nominal EHSS=20
  inf=c(4, 16)
  info=mixbeta(inf=c(1, inf))
  prior.mean=inf[1]/sum(inf)
 
# robust mixture
  mix50 <-mixbeta(informative=c(.5, inf), vague=c(.5, rob))
# emprirical Bayes power prior 
  pp <-as.powerprior(info)
## emprirical Bayes commensurate prior (takes long) 
#  cp <-as.commensurateEB(info)

 plotECSS(priorlist=list(informative=info,mixture=mix50,powerprior=pp), 
          grid=((0:40)/40), n.target=40, min.ecss=-50,
          prior.base=vague, progress="none", D=MSE)
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Error in dBetaBinomial(r, n, priormix[2, , drop = FALSE], priormix[3,  : 
#>   Assertion on 'n' failed: Must be of type 'integerish', but element 1 is not close to an integer.
#> Warning in FUN(X[[i]], ...): Uniroot() didnt coverge

#> Warning in FUN(X[[i]], ...): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in my.spg(par = 0.005, fn = lik.d, lower = 0, upper = 1, control =
#> list(maximize = TRUE, : convergence tolerance satisified at intial parameter
#> values.
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning in dbinom(y, size = k, data.mean): NaNs produced
#> Warning in uniroot(f, lower = -n.target.i + 1, upper = n.max - n.target.i, :
#> NA/Inf replaced by maximum positive value
#> Warning: Removed 41 rows containing missing values (`geom_line()`).
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# References

Wiesenfarth, M., & Calderazzo, S. (2020). Quantification of prior impact
in terms of effective current sample size. Biometrics, 76(1), 326-336.
