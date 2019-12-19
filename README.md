# SAMBA
Health research using data from electronic health records (EHR) has gained
popularity, but misclassification of EHR-derived disease status and lack of
representativeness of the study sample can result in substantial bias in
effect estimates and can impact power and type I error for association
tests. Here, the assumed target of inference is the relationship between
binary disease status and predictors modeled using a logistic regression
model. SAMBA implements several methods for obtaining bias-corrected
point estimates along with valid standard errors as proposed in Beesley and
Mukherjee (2020), currently under review.
# Installation
`SAMBA` can be downloaded from Github via the R Package `devtools`
```{r}
devtools::install_github("umich-cphds/SAMBA", build_opts = c())
```
# Vignette
Once you have `SAMBA` installed, you can type
```r
vignette("UsingSAMBA")
```
in R to bring up a tutorial on `SAMBA` and how to use it.
# Questions
For questions and comments about the implementation, please contact Alexander
Rix (alexrix@umich.edu). For questions about the method, contact Lauren Beesley
(lbeesley@umich.edu).
