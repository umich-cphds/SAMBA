#' approxdist
#'
#' This function can be used to estimate parameters in the disease model given
#' a previously-estimated marginal sensitivity. This estimation is based on
#' approximating the distribution of Dstar given Z.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population
#' of interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
#'
#' Notation:
#' \describe{
#'     \item{D}{Binary disease status of interest.}
#'     \item{Dstar}{Observed binary disease status. Potentially a misclassified
#'                  version of D. We assume D = 0 implies Dstar = 0.}
#'     \item{S}{Indicator for whether patient from population of interest is
#'              included in the analytical dataset.}
#'     \item{Z}{Covariates in disease model of interest.}
#'     \item{W}{Covariates in model for patient inclusion in analytical dataset
#'              (selection model).}
#'     \item{X}{Covariates in model for probability of observing disease given
#'              patient has disease (sensitivity model).}
#' }
#' Model Structure:
#' \describe{
#' \item{Disease Model}{\deqn{logit(P(D=1|X)) = theta_0 + theta_Z Z}}
#' \item{Selection Model}{\deqn{P(S=1|W,D)}}
#' \item{Sensitivity Model}{\deqn{logit(P(D^*=1|D=1,X)) = beta_0 + beta_X X}}
#' }
#' @param Z Matrix  with covariates in disease model
#' @param Dstar Numeric vector containing observed disease status. Should be
#'     coded as 0/1
#' @param weights Optional vector of subject-specific weights used for
#'     selection bias adjustment. Default is NULL
#' @param c_marg Marginal sensitivity, P(D^* | D)
#' @return a list with two elements: 'param', a vector with parameter estimates
#'     for disease model(intercept, logOR of Z), and 'variance', a vector of
#'     variance estimates for disease model
#' @export
approxdist <- function(Z, Dstar, c_marg, weights = NULL)
{
    if (is.data.frame(Z))
        Z <- as.matrix(Z)
    if (!is.numeric(Z))
        stop("'Z' should be a numeric matrix.")

    if (is.vector(Z))
        Z <- as.matrix(Z)
    if (!is.matrix(Z))
        stop("'Z' should be a matrix or data.frame.")

    if (!is.numeric(Dstar) || !is.vector(Dstar))
        stop("'Dstar' must be a numeric vector.")
    if (length(setdiff(0:1, unique(Dstar))) != 0)
        stop("'Dstar' must be coded 0/1.")

    if (c_marg > 1)
        stop(paste("'c_marg' greater than 1. Try estimating sensitivity",
                   "as a function of covariates."))

    n <- length(Dstar)
    if (is.null(weights)) {
        weights <- rep(1 / n, n)
    } else {
        if (length(weights) != length(Dstar))
            stop("The length of 'weights' must match the length of 'Dstar'.")
        if (!is.numeric(weights) || !is.vector(weights))
            stop("'weights' must be a numeric vector.")
        if (any(weights < 0))
            stop("'weights' must be nonnegative.")
    }

    if (length(unique(weights)) == 1) {
        fit <- stats::glm(Dstar ~ Z, family = stats::binomial())
    } else {
        # estimate using survey package due to weights
        colnames(Z) <- paste0("Z", 1:ncol(Z))
        design   <- survey::svydesign(ids = 1:n, data = data.frame(Z, Dstar),
                                      weights = weights)

        form <- paste("Dstar ~", paste(colnames(Z), collapse = " + "))
        fit  <- survey::svyglm(stats::formula(form), family = stats::binomial(),
                               design = design)
    }
    param.uc <- stats::coef(fit)[-1]
    var.uc   <- diag(summary(fit)$cov.scaled)[-1]

    p.star     <- sum(Dstar * weights) / sum(weights)
    correction <- (c_marg * (1 - p.star)) / (c_marg - p.star)

    list(param = param.uc * correction, variance = var.uc * correction ^ 2)
}
