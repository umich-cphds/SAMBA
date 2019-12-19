#' Estimate sensitivity
#'
#' \code{sensitivity} estimates (1) marginal sensitivity and (2) sensitivity as
#' a function of covariates X for a misclassified binary outcome.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates \eqn{Z} using a logistic regression model. However, \eqn{D}
#' may be misclassified, and our observed data may not well-represent the
#' population of interest. In this setting, we estimate parameters from the
#' disease model using the following modeling framework.
#'
#' Notation:
#' \describe{
#'     \item{D}{Binary disease status of interest.}
#'     \item{D*}{Observed binary disease status. Potentially a misclassified
#'                  version of D. We assume D = 0 implies D* = 0.}
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
#' \item{Sensitivity Model}{\deqn{logit(P(D* = 1| D = 1, S = 1, X)) = beta_0 + beta_X X}}
#' }
#' @param Dstar Numeric vector containing observed disease status. Should be
#'     coded as 0/1
#' @param X Numeric matrix with covariates in sensitivity model. Set to NULL
#'     to fit model with no covariates in sensitivity model. 'X' should not
#'     contain an intercept
#' @param prev marginal disease prevalence \eqn{P(D = 1)} or patient-specific
#'     \eqn{P(D = 1|X)} in population
#' @param r (optional) marginal sampling ratio, \eqn{P(S = 1|D = 1) / P(S = 1|D = 0)}.
#'     Only one of 'r' and 'weights' can be specified. Default is `NULL`
#' @param weights Optional vector of patient-specific weights used for
#'     selection bias adjustment. Only one of r and weights can be specified.
#'     Default is `NULL`
#' @return a list with two elements: (1) `c_marg`, marginal sensitivity estimate
#'     \eqn{P(D* = 1|D = 1, S = 1)}, and (2) `c_X`, sensitivity as a function of
#'     X \eqn{P(D* = 1| D = 1, S = 1, X)}
#' @export
sensitivity <- function(Dstar, X, prev, r = NULL, weights = NULL)
{
    if (!is.numeric(Dstar))
        stop("'Dstar' must be a numeric vector.")
    if (length(setdiff(0:1, unique(Dstar))) != 0)
        stop("'Dstar' must be coded 0/1.")

    n <- length(Dstar)
    if (!is.null(X)) {
        if (is.data.frame(X))
            X <- as.matrix(X)
        if (!is.numeric(X))
            stop("'X' must be numeric.")
        if (is.vector(X))
            X <- as.matrix(X)
        if (!is.matrix(X))
            stop("'X' must be a data.frame or matrix.")
        if (n != nrow(X))
            stop("The number of rows of 'X' must match the length of 'Dstar'.")
    }

    if (!is.numeric(prev) || !is.vector(prev))
        stop("'prev' must be a numeric vector.")

    if (length(prev) != 1 && length(prev) != n)
        stop("'prev' must have unit length or the same length as 'Dstar'.")

    if (length(prev) > 1)
        message('Using average prevalence to calculate marginal c')

    if (!is.null(r) && !is.null(weights))
        stop("Only one of 'r' and weights can be non-NULL.")

    check.weights(weights, n)
    if (is.null(weights))
        weights <- rep(1, n)

    if (is.null(r)) {
        r <- 1
    } else {
        if (!is.numeric(r))
            stop("'r' must be a numeric vector.")
        if (length(r) != 1 && length(r) != n)
            stop("'r' must have unit length, or the same length as 'Dstar'.")
    }

    p.star <- sum(Dstar * weights) / sum(weights)
    prevr  <- r * prev / (r * prev + 1 - prev)
    c_marg <- p.star / (r * mean(prev) / (r * mean(prev) + 1 - mean(prev)))

    c_marg <- ifelse(c_marg > 1, 1, c_marg)

    fit.beta  <- stats::glm(Dstar ~ X, family = stats::binomial(),
                            weights = weights)

    starting <- c(logit(c_marg), stats::coef(fit.beta)[-1])
    fit.sens <- try(stats::glm(Dstar ~ X, start = starting, weights = weights,
                               family = stats::binomial(modLinkprev(prevr))),
                    silent = TRUE)

    if (class(fit.sens)[1] != "try-error") {
        c1 <- expit(cbind(1, X) %*% stats::coef(fit.sens))
    } else {
        c1 <- stats::predict(fit.beta, type = 'response') / prevr
    }

    c_marg <- ifelse(c_marg < 1, c_marg, NA)
    c1     <- ifelse(c1 > 1, rep(1, length(Dstar)), c1)
    list(c_marg = c_marg, c_X = as.vector(c1))
}
