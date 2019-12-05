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
#' @param Z matrix or data frame with covariates in disease model
#' @param X matrix or data frame with covariates in sensitivity model. Set to
#'     NULL to fit model with no covariates in sensitivity model.
#' @param Dstar matrix or data frame containing observed disease status
#' @param param_current vector of starting values for theta and beta
#'     (theta, beta). Theta is the parameter of the disease model, and beta is
#'     the parameter of the sensitivity model.
#' @param beta0_fixed (optional) fixed value for the intercept of the
#'     sensitivity model.
#' @param weights (optional) vector of subject-specific weights used for
#'     selection bias adjustment.
#' @param tol stop estimation when subsequent log-likelihood estimates are
#'     within this value
#' @param maxit maximum number of iterations of the estimation algorithm
#' @return param vector with parameter estimates organized as (theta, beta)
#' @return param.seq matrix containing estimated parameter values across
#'     iterations of the expectation algorithm
#' @return loglik.seq vector of log-likelihood values across iterations of the
#'     expectation algorithm
#' @export
obsloglikEM <- function(Z, X, Dstar, param_current, beta0_fixed = NULL,
                        weights = NULL, tol = 1e-6, maxit = 50)
{
    X <- as.matrix(X)
    Z <- as.matrix(Z)

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
    if (is.null(weights))
        weights <- rep(1, nrow(X))
    while (!converged && it < maxit) {
        if (is.null(beta0_fixed)) {
            fit.beta <- stats::glm(Dstar ~ X, weights = p * weights,
                                   family = stats::binomial())
        } else {
            fit.beta <- stats::glm(Dstar ~ 0 + X, weights = p * weights,
                                   offset = rep(beta0_fixed, length(p)),
                                   family = stats::binomial())
        }
        fit.theta <- stats::glm(p ~ Z, family = stats::binomial(),
                                weights = weights)

        pred1 <- stats::predict(fit.theta, type = 'response')
        pred2 <- stats::predict(fit.beta, type = 'response')
        p     <- calculate.p(pred1, pred2)

        loglik <- sum(weights * Dstar * log(pred1 * pred2) +
                      weights * (1 - Dstar) * log(1 - (pred1 * pred2))
        )
        loglik.seq <- c(loglik.seq, loglik)

        it <- it + 1
        if (abs(loglik.seq[it] - loglik.seq[it - 1]) < tol)
            converged <- TRUE

        p <- c(stats::coef(fit.theta), stats::coef(fit.beta))
        param.seq <- rbind(param.seq,  p)
    }

    param <- c(stats::coef(fit.theta), beta0_fixed, stats::coef(fit.beta))
    list(param = param, param.seq = param.seq, loglik.seq = loglik.seq)
}
