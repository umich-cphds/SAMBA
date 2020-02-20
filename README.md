<!-- Badges -->
[![CRAN
Version](https://img.shields.io/cran/v/mianet?style=flat-square&color=blue&label=CRAN)](https://cran.r-project.org/package=SAMBA)
[![GitHub
Release](https://img.shields.io/github/v/release/umich-cphds/SAMBA?include_prereleases&label=Github&style=flat-square&color=blue)](https://github.com/umich-cphds/SAMBA)
[![Travis
CI](https://img.shields.io/travis/umich-cphds/SAMBA?style=flat-square)](https://travis-ci.org/umich-cphds/SAMBA)

SAMBA
=====

Health research using data from electronic health records (EHR) has
gained popularity, but misclassification of EHR-derived disease status
and lack of representativeness of the study sample can result in
substantial bias in effect estimates and can impact power and type I
error for association tests. Here, the assumed target of inference is
the relationship between binary disease status and predictors modeled
using a logistic regression model. SAMBA implements several methods for
obtaining bias-corrected point estimates along with valid standard
errors as proposed in Beesley and Mukherjee (2020), currently under
review.

Installation
============

`SAMBA` can be downloaded from CRAN via the R Package `devtools`

    install.packages("SAMBA")

or downloaded from Github via the R Package `devtools`

    devtools::install_github("umich-cphds/SAMBA", build_opts = c())

The Github version may contain bug fixes not yet present on CRAN, so if
you are experiencing issues, you may want to try the Github version of
the package. \# Vignette Once you have `SAMBA` installed, you can type

    vignette("UsingSAMBA")

in R to bring up a tutorial on `SAMBA` and how to use it.

Questions
=========

For questions and comments about the implementation, please contact
Alexander Rix (<alexrix@umich.edu>). For questions about the method,
contact Lauren Beesley (<lbeesley@umich.edu>).

Reference
=========

Statistical inference for association studies using electronic health
records: handling both selection bias and outcome misclassification
Lauren J Beesley, Bhramar Mukherjee medRxiv
[2019.12.26.19015859](https://doi.org/10.1101/2019.12.26.19015859)
