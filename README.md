<!-- badges: start -->
[![Github
Version](https://img.shields.io/badge/Github-0.9.0-informational.svg?style=flat)](https://github.com/umich-cphds/SAMBA)
[![Travis
CI](https://travis-ci.org/umich-cphds/SAMBA.svg?branch=master)](https://travis-ci.org/umich-cphds/SAMBA)
<!-- badges: end -->

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

`SAMBA` can be downloaded from Github via the R Package `devtools`

    devtools::install_github("umich-cphds/SAMBA", build_opts = c())

Vignette
========

Once you have `SAMBA` installed, you can type

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
