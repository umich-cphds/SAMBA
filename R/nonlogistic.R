#' nonlogistic
#' This function can be used to estimate parameters in the disease model given
#' previously-estimated sensitivity. This function fits a logistic regression
#' model for D given Z using a non-logistic link function for Dstar given Z and
#' sensitivity.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population
#' of interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
#' Notation:
#' D = binary disease status of interest
#' Dstar = observed binary disease status. Potentially a misclassified version
#' of D. We assume D = 0 implies Dstar = 0.
#' S = indicator for whether patient from population of interest is included in
#' the analytical dataset
#' Z = covariates in disease model of interest
#' W = covariates in model for patient inclusion in analytical dataset
#'     (selection model)
#' X = covariates in model for probability of observing disease given patient
#'     has disease (sensitivity model)
#' Model Structure:
#' Disease Model: logit(P(D=1|X)) = theta_0 + theta_Z Z
#' Selection Model: P(S=1|W,D)
#' Sensitivity Model: logit(P(Dstar=1|D=1,X)) = beta_0 + beta_X X
#' @param Z matrix or data frame with covariates in disease model
#' @param X matrix or data frame with covariates in sensitivity model. Set to
#'     NULL to fit model with no covariates in sensitivity model.
#' @param Dstar matrix or data frame containing observed disease status
#' @param weights (optional) vector of subject-specific weights used for
#'     selection bias adjustment.
#' @param c_X sensitivity as a function of X, P(observe disease | have disease, X)
#' @return param vector with parameter estimates for disease model
#'   (intercept, logOR of Z)
#' @return variance vector of variance estimates for disease model
#' @export
nonlogistic <- function(Z, X, Dstar, weights = NULL, c_X)
{
    Z = data.frame(Z)
    c_X  <- as.vector(c_X)

    fitTheta <- glm(Dstar ~ as.matrix(Z), family = binomial())
    starting <- as.numeric(coef(fitTheta))
    link <- modLink(c_X)
    if (is.null(weights)) {
        fit <- stats::glm(Dstar ~ as.matrix(Z), family = stats::binomial(link = link), start = starting)

        var <- diag(summary(fit)$cov.scaled)
    } else {
        ### Estimation using survey package:
        names(Z) = paste0("Z", 1:ncol(Z))
        starting = data.frame(matrix(rep(coef(fitTheta), length(Dstar)), byrow = T, ncol = length(coef(fitTheta))))
        names(starting) = paste0('start', c(0:dim(as.matrix(Z))[2]))
        design = survey::svydesign(ids = c(1:length(Dstar)),strata = NULL, data = data.frame(Z, Dstar = Dstar, c_X = c_X, starting),  weights = weights)
        formula = paste( names(Z),collapse = '+')
        formula_starting = paste0('c(',paste( 'data$',names(starting),collapse = '[1],', sep = ''),'[1])')
        formula_long = paste0('survey::svyglm(Dstar~',formula,', family = binomial(link=link), start = ',formula_starting,', design = design)')
        fit = eval(parse(text = formula_long))

        var <- diag(vcov(fit))
    }

    param <- stats::coef(fit)
    list(param = param, variance = var)
}
