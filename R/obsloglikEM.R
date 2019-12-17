#' obsloglikEM
#' This function can be used to jointly estimate disease model and sensitivity
#' model parameters using an expectation-maximization algorithm.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population
#' of interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
#' Notation:
#' D = binary disease status of interest
#' Dstar = observed binary disease status. Potentially a misclassified version
#' of D. We assume D=0 implies Dstar=0.
#' S = indicator for whether patient from population of interest is included
#' in the analytical dataset
#' Z = covariates in disease model of interest
#' W = covariates in model for patient inclusion in analytical dataset
#' (selection model)
#' X = covariates in model for probability of observing disease given patient
#' has disease (sensitivity model)
#' Model Structure:
#' Disease Model: logit(P(D=1|X)) = theta_0 + theta_Z Z
#' Selection Model: P(S=1|W,D)
#' Sensitivity Model: logit(P(Dstar=1|D=1,X)) = beta_0 + beta_X X
#'
#' @param Dstar matrix or data frame containing observed disease status
#' @param Z matrix or data frame with covariates in disease model
#' @param X Numeric matrix with covariates in sensitivity model. Set to NULL
#'     to fit model with no covariates in sensitivity model. 'X' should not
#'     contain an intercept.
#' @param param_current vector of starting values for theta and beta
#'     (theta, beta). Theta is the parameter of the disease model, and beta is
#'     the parameter of the sensitivity model.
#' @param beta0_fixed (optional) fixed value for the intercept of the
#'     sensitivity model.
#' @param weights (optional) vector of subject-specific weights used for
#'     selection bias adjustment.
#' @param expected Whether or not to calculate the covariance matrix via the
#'     expected fisher information matrix. Default is TRUE
#' @param tol stop estimation when subsequent log-likelihood estimates are
#'     within this value
#' @param maxit maximum number of iterations of the estimation algorithm
#' @return param vector with parameter estimates organized as (theta, beta)
#' @return param.seq matrix containing estimated parameter values across
#'     iterations of the expectation algorithm
#' @return loglik.seq vector of log-likelihood values across iterations of the
#'     expectation algorithm
#' @export
obsloglikEM <- function(Dstar, Z, X, param_current, beta0_fixed = NULL,
                        weights = NULL, expected = TRUE, tol = 1e-6, maxit = 50)
{
    if (is.data.frame(Z))
        Z <- as.matrix(Z)
    if (!is.numeric(Z))
        stop("'Z' should be a numeric matrix.")

    if (is.vector(Z))
        Z <- as.matrix(Z)
    if (!is.matrix(Z))
        stop("'Z' should be a matrix or data.frame.")

    if (!is.null(X)) {
        if (!is.numeric(X))
            stop("'X' must be numeric.")
        if (is.vector(X))
            X <- as.matrix(X)
        if (!is.matrix(X))
            stop("'X' must be a data.frame or matrix.")
    }

    if (!is.numeric(Dstar) || !is.vector(Dstar))
        stop("'Dstar' must be a numeric vector.")
    if (length(setdiff(0:1, unique(Dstar))) != 0)
        stop("'Dstar' must be coded 0/1.")

    n <- length(Dstar)
    if (nrow(Z) != n)
        stop("The number of rows of 'Z' must match the length of 'Dstar'.")
    if (!is.null(X) && nrow(X) != n)
        stop("The number of rows of 'X' must match the length of 'Dstar'.")

    check.weights(weights, n)
    if (is.null(weights))
        w <- rep(1, n)
    else
        w <- weights

    if (!is.logical(expected) || length(expected) > 1)
        stop("'expected' must be a length one logical.")

    # initialise p for EM
    theta <- param_current[1:(1 + ncol(Z))]
    beta  <- param_current[-(1:(1 + ncol(Z)))]
    pred1 <- expit(cbind(1, Z) %*% theta)
    pred2 <- expit(cbind(1, X) %*% beta)

    calculate.p <- function(pred1, pred2)
    {
        Dstar + (1 - Dstar) * ((pred1 * (1 - pred2)) /
                              (pred1 * (1 - pred2) + (1 - pred1)))
    }
    p <- calculate.p(pred1, pred2)

    it <- 1
    converged <- F

    param.seq  <- c(theta, beta0_fixed, beta)
    loglik.seq <- -10 ^ 9
    while (!converged && it < maxit) {
        if (is.null(beta0_fixed)) {
            fit.beta <- stats::glm(Dstar ~ X, weights = p * w,
                                   family = stats::binomial())
        } else {
            fit.beta <- stats::glm(Dstar ~ 0 + X, weights = p * w,
                                   offset = rep(beta0_fixed, length(p)),
                                   family = stats::binomial())
        }
        fit.theta <- stats::glm(p ~ Z, family = stats::binomial(),
                                weights = weights)

        pred1 <- stats::predict(fit.theta, type = 'response')
        pred2 <- stats::predict(fit.beta, type = 'response')
        p     <- calculate.p(pred1, pred2)

        loglik <- sum(w * Dstar * log(pred1 * pred2) +
                      w * (1 - Dstar) * log(1 - (pred1 * pred2))
        )
        loglik.seq <- c(loglik.seq, loglik)

        it <- it + 1
        if (abs(loglik.seq[it] - loglik.seq[it - 1]) < tol)
            converged <- TRUE

        p <- c(stats::coef(fit.theta), stats::coef(fit.beta))
        param.seq <- rbind(param.seq,  p)
    }

    param <- c(stats::coef(fit.theta), beta0_fixed, stats::coef(fit.beta))
    theta <- param[1:(ncol(Z) + 1)]
    beta  <- param[-(1:(ncol(Z) + 1))]

    if (is.null(weights)) {
        var <- obsloglik_var(Dstar, Z, X, theta, beta, beta0_fixed, expected)
    } else {
        var <- obsloglik_var_weighted(Dstar, Z, X, theta, beta, beta0_fixed,
                                      weights, expected)
    }

    structure(list(param = param, var = var, param.seq = param.seq, loglik.seq =
                   loglik.seq, Dstar = Dstar, X = X, Z = Z, weights = weights,
                   beta0_fixed = beta0_fixed), class = "SAMBA.fit")

}
