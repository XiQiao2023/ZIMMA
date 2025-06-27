
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Decoding Microbiome Dual Mediation: Introducing ZIMMA for Enhanced Zero-Inflated Data Analysis

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

This example demonstrates how the dual microbiome mediation method can
be applied to a longitudinal study investigating gut microbiome
structure following hematopoietic stem cell transplantation (HCT)
([Schluter et al.,
2020](https://www.nature.com/articles/s41586-020-2971-8)). The
processced 16S data is public available through [Liao et
al.,2021](https://www.nature.com/articles/s41597-021-00860-8).

We analyzed data from a cohort of $N = 256$ patients who underwent a
single HCT and had microbiome profiles available on the day of
transplantation (day 0). We included amplicon sequence variants (ASVs)
with a prevalence of at least 10%, resulting in $p = 1299$ ASVs spanning
21 classes, 34 orders, 66 families, and 231 genera across 13 phyla.

``` r
library(ZIMMA)

data("Hema_physeq")
Hema_physeq
#> Loading required package: phyloseq
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 3911 taxa and 256 samples ]
#> sample_data() Sample Data:       [ 256 samples by 36 sample variables ]
#> tax_table()   Taxonomy Table:    [ 3911 taxa by 6 taxonomic ranks ]

## Appling taxa filtering procedure: Prevalence > 1%
M <- as.data.frame(t(Hema_physeq@otu_table))
M <- M[, colMeans(M != 0) > 0.01] 

## Count the number of unique features at each taxonomic level
taxa_to_keep <- colnames(M)
physeq_filtered <- prune_taxa(taxa_to_keep, Hema_physeq)
taxa_table_filtered <- tax_table(physeq_filtered)
```

## Apply the dual mediation framework using ZIMMA()

We hypothesized that the gut microbiome mediates the association between
piperacillin/tazobactam administration prior to HCT and post-HCT white
blood cell counts in transplantation patients. Confounding variables for
the mediator model included hematopoietic cell source, underlying
disease, and other medications administered prior to HCT. In addition to
these, the outcome model also adjusted for medications administered
after HCT.

``` r
Hema_Phylum_Result = ZIMMA(object = physeq_filtered, 
               Outcome = "log_BC", 
               Treat = "Drug_piperacillinOrtazobactam_Prior_Use", 
               C.med = c("HCTSource","Disease",
                         "Drug_ciprofloxacin_Prior_Use",
                         "Drug_ganciclovir_Prior_Use"), 
               C.out =  c("HCTSource","Disease",
                          "Drug_ciprofloxacin_Prior_Use",
                          "Drug_ganciclovir_Prior_Use",
                          "Drug_imipenem_cilastatin_Post_Use",
                          "Drug_piperacillinOrtazobactam_Post_Use",
                          "Drug_atovaquone_Post_Use",
                          "Drug_isavuconazonium_sulfate_Post_Use",
                          "Drug_meropenem_Post_Use"), 
               M.level = "Phylum",
               n_iter = 20000, burn_in = 5000,
               Size = "RLEpseudo", pseudo = 0.5)

save(Hema_Phylum_Result,file = "Hema_Phylum_Result.rda")
```

## Get the summary of the posterior samples using ZIMMA.summary()

``` r
sum = ZIMMA.summary(Hema_Phylum_Result)
```

## Select the active mediator carrying the indirect effect through abundance

Active Threshold: PIP \> 0.5

``` r
NIE.A = which(sum$betaT[,1]>0.5 & sum$alphaM[,1]>0.5) 
NIE.A
#> <not present> 
#>             3
# get PIP, posterior mean, posterior median, HPD
sum$betaT[NIE.A,, drop = FALSE]
#>                     PIP      mean    median HPD.Lower HPD.Upper
#> <not present> 0.9548667 -2.175247 -2.172305 -3.354904 -1.178494
sum$alphaM[NIE.A,,drop = FALSE]
#>                     PIP      Mean    Median   HPD.Lower HPD.Upper
#> <not present> 0.5874667 0.1178958 0.1155591 -0.02453327 0.3099058
```

## Select the active mediator carrying the indirect effect through prevalence

Active Threshold: PIP \> 0.5

``` r
NIE.P = which(sum$gammaT[,1]>0.5 & sum$alphaM[,1]>0.5)
NIE.P
#> named integer(0)
sum$gammaT[NIE.P,,drop = FALSE]
#>      PIP mean median HPD.Lower HPD.Upper
sum$alphaM[NIE.P,,drop = FALSE]
#>      PIP Mean Median HPD.Lower HPD.Upper
```
