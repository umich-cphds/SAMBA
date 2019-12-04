#' approxdist
#'
#' This function can be used to estimate parameters in the disease model given
#' previously-estimated marginal sensitivity. This estimation is based on
#' approximating the distribution of Dstar given Z.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population
#' of interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
#'
#' Notation:
#' D = binary disease status of interest
#' Dstar = observed binary disease status. Potentially a misclassified version
#' of D. We assume D = 0 implies Dstar = 0. S = indicator for whether patient
#' from population of interest is included in the analytical dataset
#' Z = covariates in disease model of interest
#' W = covariates in model for patient inclusion in analytical dataset
#' (selection model)
#' X = covariates in model for probability of observing disease given patient
#' has disease (sensitivity model)
#' Model Structure:
#' Disease Model: logit(P(D=1|X)) = theta_0 + theta_Z Z
#' Selection Model: P(S=1|W,D)
#' Sensitivity Model: logit(P(Dstar=1|D=1,X)) = beta_0 + beta_X X
#' @param Z Matrix or data frame with covariates in disease model
#' @param Dstar Matrix or data frame containing observed disease status
#' @param weights Optional vector of subject-specific weights used for
#'     selection bias adjustment. Default is NULL
#' @param c_marg Marginal sensitivity, P( observe disease | have disease )
#' @return a list with two elements: param, a vector with parameter estimates
#'     for disease model(intercept, logOR of Z), and variance, a vector of
#'     variance estimates for disease model
#' @export
approxdist <- function(Z, Dstar, c_marg, weights = NULL)
{
    if (is.matrix(Z))
        Z <- as.data.frame(Z)
    if (!is.data.frame(Z))
        stop("'Z' should be a matrix or data.frame.")
    if (!is.numeric(Dstar))
        stop("'Dstar' should be a binary variable.")

    if (c_marg > 1)
        stop(paste("'c_marg greater than 1. Try estimating sensitivity",
                   "as a function of covariates."))
    n <- length(Dstar)
    if (is.null(weights))
        weights <- rep(1 / n, n)

    if (length(unique(weights)) == 1) {
        fit <- stats::glm(Dstar ~ as.matrix(Z), family = stats::binomial())
    } else {
        # estimate using survey package due to weights
        names(Z) <- paste0("Z", 1:ncol(Z))
        design   <- survey::svydesign(ids = 1:n, data = data.frame(Z, Dstar),
                                      weights = weights)

        f   <- formula(paste("Dstar ~", paste(names(Z), collapse = " + ")))
        fit <- survey::svyglm(f, family = stats::binomial(), design = design)
    }
    param.uc <- stats::coef(fit)[-1]
    var.uc   <- diag(summary(fit)$cov.scaled)[-1]

    p.star     <- sum(Dstar * weights) / sum(weights)
    correction <- (c_marg * (1 - p.star)) / (c_marg - p.star)

    list(param = param.uc * correction, variance = var.uc * correction ^ 2)
}
