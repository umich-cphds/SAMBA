#' Estimate parameters in the disease model given a previously-estimated
#' sensitivity.
#'
#' \code{nonlogistic} fits a logistic regression model for D given Z using a
#' non-logistic link function for Dstar given Z and sensitivity.
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
#' @param Z matrix or data frame with covariates in disease model
#' @param Dstar matrix or data frame containing observed disease status
#' @param weights (optional) vector of subject-specific weights used for
#'     selection bias adjustment.
#' @param c_X sensitivity as a function of X, P(observe disease | have disease, X)
#' @return param vector with parameter estimates for disease model
#'   (intercept, logOR of Z)
#' @return variance vector of variance estimates for disease model
#' @export
nonlogistic <- function(Z, Dstar, c_X, weights = NULL)
{
    if (is.data.frame(Z))
        Z <- as.matrix(Z)
    if (!is.numeric(Z))
        stop("'Z' must be a numeric matrix.")
    if (is.vector(Z))
        Z <- as.matrix(Z)
    if (!is.matrix(Z))
        stop("'Z' must be a numeric matrix.")

    if (!is.numeric(Dstar) || !is.vector(Dstar))
        stop("'Dstar' must be a numeric vector.")
    if (length(setdiff(0:1, unique(Dstar))) != 0)
        stop("'Dstar' must be coded 0/1.")

    if (nrow(Z) != length(Dstar))
        stop("The number of rows of 'Z' must match the length of 'Dstar'.")

    if (!is.numeric(c_X) || !is.vector(c_X))
        stop("'c_X' must be a numeric vector.")
    if (length(c_X) != length(Dstar))
        stop("'c_X' must have the same length as 'Dstar'.")

    if (!is.null(weights)) {
        if (length(weights) != length(Dstar))
            stop("The length of 'weights' must match the length of 'Dstar'.")
        if (!is.numeric(weights) || !is.vector(weights))
            stop("'weights' must be a numeric vector.")
        if (any(weights < 0))
            stop("'weights' must be nonnegative.")
    }

    fitTheta <- stats::glm(Dstar ~ Z, family = stats::binomial())
    starting <- stats::coef(fitTheta)
    if (is.null(weights)) {
        fit <- stats::glm(Dstar ~ Z, start = starting, family =
                          stats::binomial(modLink(c_X)))

        var <- diag(summary(fit)$cov.scaled)
    } else {
        # Estimate using survey package if weighted
        colnames(Z) <- paste0("Z", 1:ncol(Z))
        design <- survey::svydesign(ids = 1:length(Dstar), strata = NULL,
                                    data = data.frame(Z, Dstar = Dstar),
                                    weights = weights)

        form <- paste("Dstar ~", paste(colnames(Z), collapse = "+"))
        fit  <- survey::svyglm(stats::formula(form), design = design,
                               family = stats::binomial(modLink(c_X)),
                               start = stats::coef(fitTheta))

        var <- diag(stats::vcov(fit))
    }

    param <- stats::coef(fit)
    list(param = param, variance = var)
}
