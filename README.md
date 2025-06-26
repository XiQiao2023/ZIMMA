
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ZIMMA

<!-- badges: start -->
<!-- badges: end -->

Given the complex interactions between the microbiome, the host, and
external factors, causal mediation analysis is essential for unraveling
how dysbiosis or microbial imbalance mediates the effects of
interventions or environmental exposures on health outcomes. However,
zero inflation in microbiome count data complicates high-dimensional
mediation analysis, as frequently employed zero-inflated models often
struggle to distinguish true zero inflation from the underlying count
distribution. To address this, we employ a zero-inflated negative
binomial distribution for each mediator within a Bayesian framework,
incorporating spike-and-slab priors to enable sparsity in estimating
natural indirect effects (NIE) and an informative prior based on nonzero
counts to improve dispersion estimation and control zero sources.
Recognizing the distinct biological mechanisms underlying microbial
presence versus abundance, we developed a dual mediation model, ZIMMA,
to separate the NIE into pathways for mediator abundance and prevalence.
Extensive simulations demonstrate ZIMMA’s superior performance in
capturing distinct mediation mechanisms for both rare and abundant
species compared to existing methods. Its application to human
microbiome studies examining the effects of dietary intake and metabolic
syndrome underscores its efficacy in identifying key microbial
mediators, offering valuable insights and biological interpretation into
the role of the microbiome in disease physiology and health sciences.

## Installation

You can install the development version of ZIMMA like so:

``` r
devtools::install_github("XiQiao2023/ZIMMA")
```

## Load the Phyloseq Object

This is a basic example which shows you how to solve a common problem:

``` r
library(ZIMMA)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
