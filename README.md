
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mutalisk

<!-- badges: start -->

[![R-CMD-check](https://github.com/selkamand/mutalisk/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/mutalisk/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/selkamand/mutalisk/branch/main/graph/badge.svg)](https://app.codecov.io/gh/selkamand/mutalisk?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/mutalisk)](https://CRAN.R-project.org/package=mutalisk)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

mutalisk is an R package that produces cohort-level visualisations from
the files output by [**mutalisk**](http://mutalisk.org/).

Note [CRUX](https://ccicb.shinyapps.io/crux/) is a GUI tool that makes
it easy to export data in mutalisk compatible formats, and also provides
a GUI interface to this package in the ‘Mutational Signatures’ module.
For a no-code solution, consider using the web app.

## Installation

You can install the development version of mutalisk from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("selkamand/mutalisk")
```

## Usage

First you have to generate mutalisk file:

1.  Start by following the instructions on the
    [**mutalisk**](http://mutalisk.org/) website

2.  Once the results are ready, Select
    `Mutational Signature (Best only)`

3.  Then click `Get the selected result for all samples at once`

4.  Unzip the file. This mutalisk directory can

``` r
library(mutalisk)

# Read files in mutalisk_dir into 
mutalisk_dir <- system.file("lusc_tcga",package = "mutalisk")
mutalisk_df <- mutalisk_best_signature_directory_to_dataframe(mutalisk_dir)

# Plot Cohort Data
plot_stacked_bar(mutalisk_df)
#> Lumping together signatures with contributions < 0.1 in all samples as 'Other'
#> Metadata columns: SignaturesAllContributionsBelowMin
```

<img src="man/figures/README-example-1.png" width="100%" />
