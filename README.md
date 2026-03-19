AUGglmmTMB
================

# AUGglmmTMB

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

- `mpl_fitter()`: fits penalized binomial mixed models using iterative
  algorithms based on Košuta et al. The function supports penalties on
  both fixed and random effects and provides stable estimation in the
  presence of separation and boundary estimates.

- `get_psi()`: estimates penalty parameters for the random-effects
  covariance structure using an iterative procedure consistent with the
  penalized likelihood framework.

- `get_data_plot_cloglik()`: evaluates the conditional log-likelihood
  over a grid of shrinkage parameters for the random-effects covariance
  matrix, enabling data-driven selection of the penalty strength.

## Installation

You can install the development version from GitHub:

\`\`\`r \# install.packages(“remotes”)
remotes::install_github(“rokblagus/AUGglmmTMB”)
