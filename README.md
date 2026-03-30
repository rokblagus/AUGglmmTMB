AUGglmmTMB
================

## Overview

`AUGglmmTMB` is an R package for the analysis of binomial mixed models
in settings with sparse and clustered data. It provides methods for
stable estimation in the presence of separation and boundary estimates
of random effects covariance matrices.

The package implements penalized likelihood approaches based on
extensions of Firth-type penalty for the fixed effects and Inverse
Wishart type penalty for the random effects covariance matrix through
data augmentation techniques as presented in Košuta et al. These methods
ensure finite and interpretable parameter estimates while preserving
desirable statistical properties such as parameterization invariance.

## Motivation

In many applied settings, particularly in clinical and epidemiological
research, data are sparse and structured hierarchically. Standard
maximum likelihood estimation in binomial mixed models may fail due to:

- separation in binary or count outcomes,
- rare events or rare exposures,
- boundary estimates of variance components.

These issues lead to unstable estimation and invalid inference.
`AUGglmmTMB` provides practical tools to address these challenges within
a unified framework.

## Features

- Penalized likelihood estimation for GLMMs based on Firth-type
  corrections  
- Data augmentation approaches for stable estimation  
- Handling of separation in clustered data  
- Regularization of random effects covariance parameters  
- Compatibility with `glmmTMB` workflows

## Key Functions

- `AUGglmmTMB()`: fits penalized binomial mixed models using iterative
  algorithms based on Košuta et al. The function supports penalties on
  both fixed and random effects and provides stable estimation in the
  presence of separation and boundary estimates.

- `mpl_fitter()`: fits penalized binomial mixed models using iterative
  algorithms based on Košuta et al. The function supports penalties on
  both fixed and random effects and provides stable estimation in the
  presence of separation and boundary estimates.

- `get_psi()`: estimates penalty parameters for the random-effects
  covariance structure using an iterative procedure consistent with the
  penalized likelihood framework.

- `get_data_plot_cloglik()`: evaluates the conditional log-likelihood
  over a grid of shrinkage parameters for the random-effects covariance
  matrix, supporting data-driven selection of the penalty strength
  obtained using `get_psi()`.

## Installation

You can install the development version from GitHub:

`install.packages("remotes");remotes::install_github("rokblagus/AUGglmmTMB")`

## Examples

The functionality of the package is illustrated using a dataset on
European passerine birds, originally published by Bandelj et al. The
data used in their study included additional sampling to address the
separation issue. In contrast, the dataset `data(birds)` provided with
this R package corresponds to the original dataset, prior to this
additional sampling.

First, we fit the model using ML, as implemented in the R package
`glmmTMB`, using the default settings.

``` r
library(glmmTMB)
```

    ## Warning: package 'glmmTMB' was built under R version 4.4.3

``` r
library(AUGglmmTMB)
data(birds)
fit_ml<-glmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=birds,family=binomial(link = "logit"))
```

    ## Warning in finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old): Model convergence
    ## problem; singular convergence (7). See vignette('troubleshooting'),
    ## help('diagnose')

``` r
summary(fit_ml)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## parasites ~ migration + food + (migration | phylogenetic) + (1 |      species)
    ## Data: birds
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     315.6     350.7    -148.8     297.6       357 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups       Name        Variance Std.Dev. Corr  
    ##  phylogenetic (Intercept) 0.03439  0.1854         
    ##               migration   0.75955  0.8715   -1.00 
    ##  species      (Intercept) 0.39084  0.6252         
    ## Number of obs: 366, groups:  phylogenetic, 8; species, 42
    ## 
    ## Conditional model:
    ##               Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)   -1.84896        NaN     NaN      NaN
    ## migration     -1.87632        NaN     NaN      NaN
    ## foods       -326.00895        NaN     NaN      NaN
    ## foodv          2.36462        NaN     NaN      NaN
    ## foodz         -0.07317        NaN     NaN      NaN

Due to the separation issue, the coefficient for `foods` is
unreasonable. The estimated random effects covariance matrix at the
level of `phylogenetic` is on the boundary of the parameter space, as
evident by the correlation parameter `Corr` estimated at -1.
Furtheremore, the standar errors cannot be computed for any of the
model’s parameters.

``` r
fit_ml$sdr
```

    ## sdreport(.) result
    ##            Estimate Std. Error
    ## beta    -1.84896183        NaN
    ## beta    -1.87632269        NaN
    ## beta  -326.00895307        NaN
    ## beta     2.36462228        NaN
    ## beta    -0.07317489        NaN
    ## theta   -1.68505042        NaN
    ## theta   -0.13751309        NaN
    ## theta -371.21301045        NaN
    ## theta   -0.46972634        NaN
    ## Maximum gradient component: 6.988847e-05

To solve these issues, we fit the model using MPL, using the default
penalty for the fixed effects penalty, inverse-Wishart penalty with
$\nu=3$ and $\Psi=I_2$ ($2\times 2$ identity matrix) for the random
effects covariance matrix at the level of `phylogenetic` and no penalty
for the random effects covariance by `species`. We use Algorithm 1 from
Košuta et al.

``` r
fit_mpl_1<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=birds,link = "logit",
                penOpt = AUGglmmTMBPenalty(autrepen =FALSE,nu=list(3),psi=list(diag(1,2,2))),
                      control=AUGglmmTMBControl(fit_pGLM = FALSE,maxiter = 50,tol=1e-5,save_coef=TRUE)   )
summary(fit_mpl_1$fit$fit)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## cbind(Y, M - Y) ~ -1 + X + (-1 + Z1 | grouping1) + (-1 + Z2 |      grouping2)
    ## Data: xdfa
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##        NA        NA        NA        NA       369 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups    Name          Variance Std.Dev. Corr  
    ##  grouping1 Z1(Intercept) 0.1521   0.3900         
    ##            Z1migration   0.1667   0.4083   -0.01 
    ##  grouping2 Z2            0.4047   0.6362         
    ## Number of obs: 378, groups:  grouping1, 14; grouping2, 48
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## X(Intercept) -2.03708    0.41458  -4.914 8.94e-07 ***
    ## Xmigration   -1.38678    1.12539  -1.232   0.2179    
    ## Xfoods       -1.16077    3.14559  -0.369   0.7121    
    ## Xfoodv        2.30047    1.37073   1.678   0.0933 .  
    ## Xfoodz       -0.01401    0.51790  -0.027   0.9784    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The coefficient for `foods` is now more reasonable, with the random
effects covariance matrix at the level of `phylogenetic` estimated away
from the boundary. Unlike for ML, the standard errors for model’s
parameters can now be computed. However, there is still an issue with
the SE for one of the random effects parameters, implying that the
chosen prior for the random effects is not optimal.

``` r
fit_mpl_1$fit$fit$sdr
```

    ## sdreport(.) result

    ## Warning in sqrt(diag(object$cov.fixed)): NaNs produced

    ##           Estimate Std. Error
    ## beta  -2.037080641  0.4145798
    ## beta  -1.386780220  1.1253945
    ## beta  -1.160766473  3.1455864
    ## beta   2.300469992  1.3707286
    ## beta  -0.014005846  0.5179004
    ## theta -0.941613825  0.8747518
    ## theta -0.895685689  0.9383581
    ## theta -0.008284268        NaN
    ## theta -0.452297046  0.4642174
    ## Warning:
    ## Hessian of fixed effects was not positive definite.
    ## Maximum gradient component: 0.7870051

Therefore, we use the data-driven procedure proposed by Košuta et al.,
to determine the penalty parameters for the random effects penalty.

``` r
fit_mpl_2<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=birds,link = "logit",
                penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=TRUE,ntaus=20),
                      control=AUGglmmTMBControl(fit_pGLM = FALSE,maxiter = 50,tol=1e-5,save_coef=TRUE)   )
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The value of $\nu$ used by this data-driven procedure is equal to
$\nu=2q-1$, as suggested by Košuta et al. The plot produced by
`plot=TRUE,ntaus=20`, shows the conditional likelihood as a function of
$\tau$, with the blue vertical line denoting the chosen value of $\tau$.
Using this data-driven penalty

``` r
fit_mpl_2$optre
```

    ## $tau
    ## [1] 0.126512
    ## 
    ## $psi
    ##            [,1]       [,2]
    ## [1,]  0.9291683 -0.8215102
    ## [2,] -0.8215102  4.7267192

the standard errors are now available for all model’s parameters

``` r
fit_mpl_2$fit$fit$sdr
```

    ## sdreport(.) result
    ##          Estimate Std. Error
    ## beta  -1.95674377  0.6069806
    ## beta  -1.54866333  1.4033417
    ## beta  -0.88192275  3.3835891
    ## beta   2.21984059  1.5328790
    ## beta  -0.06668198  0.5690398
    ## theta -0.97543819  0.8599558
    ## theta -0.12680044  0.9642099
    ## theta -0.44180381  3.8750655
    ## theta -0.48510956  0.4892424
    ## Maximum gradient component: 0.9157068

``` r
summary(fit_mpl_2$fit$fit)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## cbind(Y, M - Y) ~ -1 + X + (-1 + Z1 | grouping1) + (-1 + Z2 |      grouping2)
    ## Data: xdfa
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     334.3     369.8    -158.2     316.3       369 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups    Name          Variance Std.Dev. Corr  
    ##  grouping1 Z1(Intercept) 0.1421   0.3770         
    ##            Z1migration   0.7760   0.8809   -0.40 
    ##  grouping2 Z2            0.3790   0.6156         
    ## Number of obs: 378, groups:  grouping1, 14; grouping2, 48
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)   
    ## X(Intercept) -1.95674    0.60698  -3.224  0.00127 **
    ## Xmigration   -1.54866    1.40334  -1.104  0.26979   
    ## Xfoods       -0.88192    3.38359  -0.261  0.79436   
    ## Xfoodv        2.21984    1.53288   1.448  0.14757   
    ## Xfoodz       -0.06668    0.56904  -0.117  0.90671   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

We now use a computationally more efficient approximation of MPL, based
on a single iteration of Algorithm 2 in Košuta et al.

``` r
fit_mpl_3<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=birds,link = "logit",
                penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=TRUE,ntaus=20),
                      control=AUGglmmTMBControl(fit_pGLM = TRUE,maxiter = 1,tol=1e-5,save_coef=TRUE)   )
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

which gives very similar results than the exact implementation based on
Algorithm 1: the difference between the two procedures relative to the
SE obtained by the exact procedure is smaller than 1% of the estimated
SE for all model’s parameters

``` r
(fit_mpl_2$fit$fit$sdr$par.fixed-fit_mpl_3$fit$fit$sdr$par.fixed)/sqrt(diag(fit_mpl_2$fit$fit$sdr$cov.fixed))*100
```

    ##         beta         beta         beta         beta         beta        theta 
    ## -0.033364739  0.027908284  0.006830284 -0.017427823 -0.030092842  0.140516868 
    ##        theta        theta        theta 
    ##  0.030918940  0.017866347 -0.004162099

``` r
fit_mpl_3$fit$fit$sdr
```

    ## sdreport(.) result
    ##          Estimate Std. Error
    ## beta  -1.95654126  0.6057950
    ## beta  -1.54905498  1.4020412
    ## beta  -0.88215386  3.3831690
    ## beta   2.22010774  1.5326498
    ## beta  -0.06651074  0.5686943
    ## theta -0.97664658  0.8605040
    ## theta -0.12709856  0.9624175
    ## theta -0.44249614  3.8629408
    ## theta -0.48508920  0.4891596
    ## Maximum gradient component: 0.9150307

``` r
summary(fit_mpl_3$fit$fit)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## cbind(Y, M - Y) ~ -1 + X + (-1 + Z1 | grouping1) + (-1 + Z2 |      grouping2)
    ## Data: xdfa
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     334.3     369.8    -158.2     316.3       369 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups    Name          Variance Std.Dev. Corr  
    ##  grouping1 Z1(Intercept) 0.1418   0.3766         
    ##            Z1migration   0.7755   0.8806   -0.40 
    ##  grouping2 Z2            0.3790   0.6156         
    ## Number of obs: 378, groups:  grouping1, 14; grouping2, 48
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)   
    ## X(Intercept) -1.95654    0.60580  -3.230  0.00124 **
    ## Xmigration   -1.54905    1.40204  -1.105  0.26922   
    ## Xfoods       -0.88215    3.38317  -0.261  0.79429   
    ## Xfoodv        2.22011    1.53265   1.449  0.14747   
    ## Xfoodz       -0.06651    0.56869  -0.117  0.90690   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

To implement a penalty on both random effects, we can e.g. use

``` r
fit_mpl_4<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=birds,link = "logit",
                penOpt = AUGglmmTMBPenalty(autrepen =FALSE,nu=list(3,1),psi=list(fit_mpl_3$optre$psi,matrix(1,1,1))),
                      control=AUGglmmTMBControl(fit_pGLM = TRUE,maxiter = 1,tol=1e-5,save_coef=TRUE)   )
```

where we used the data-driven penalty for the covariance matrix at the
level of `phylogenetic` and the inverse-Wishart prior with $\nu=1$ and
$\Psi=1$ for the covariance at the level of `species`.

``` r
summary(fit_mpl_4$fit$fit)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## cbind(Y, M - Y) ~ -1 + X + (-1 + Z1 | grouping1) + (-1 + Z2 |      grouping2)
    ## Data: xdfa
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     334.4     369.8    -158.2     316.4       369 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups    Name          Variance Std.Dev. Corr  
    ##  grouping1 Z1(Intercept) 0.1419   0.3767         
    ##            Z1migration   0.7768   0.8814   -0.41 
    ##  grouping2 Z2            0.3523   0.5935         
    ## Number of obs: 378, groups:  grouping1, 14; grouping2, 48
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)   
    ## X(Intercept) -1.94664    0.60921  -3.195   0.0014 **
    ## Xmigration   -1.54581    1.41090  -1.096   0.2732   
    ## Xfoods       -0.87966    3.37982  -0.260   0.7947   
    ## Xfoodv        2.20928    1.52384   1.450   0.1471   
    ## Xfoodz       -0.06802    0.56237  -0.121   0.9037   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
fit_mpl_4$fit$fit$sdr
```

    ## sdreport(.) result
    ##          Estimate Std. Error
    ## beta  -1.94664034  0.6092142
    ## beta  -1.54580765  1.4108966
    ## beta  -0.87966236  3.3798197
    ## beta   2.20927865  1.5238398
    ## beta  -0.06802475  0.5623721
    ## theta -0.97623609  0.8563728
    ## theta -0.12627366  0.9765183
    ## theta -0.44361494  3.9169984
    ## theta -0.52166366  0.5027700
    ## Maximum gradient component: 0.9081954
