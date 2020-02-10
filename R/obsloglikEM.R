#' Estimate parameters in the disease model using observed data log-likelihood
#' using the expectation-maximization algorithm
#'
#' \code{obsloglikEM} jointly estimates the disease model and sensitivity
#' model parameters using profile likelihood methods. Estimation involves
#' an expectation-maximization algorithm.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population
#' of interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
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
#' \item{Sensitivity Model}{\deqn{logit(P(D* = 1| D = 1, S = 1, X)) = beta_0 +
#'                                    beta_X X}}
#' }
#' @param Dstar Numeric vector containing observed disease status. Should be
#'     coded as 0/1
#' @param Z Numeric matrix of covariates in disease model. 'Z' should not
#'     contain an intercept
#' @param X Numeric matrix of covariates in sensitivity model. Set to
#'     NULL to fit model with no covariates in sensitivity model. 'X' should not
#'     contain an intercept
#' @param start Numeric vector of starting values for theta and beta
#'     (theta, beta). Theta is the parameter of the disease model, and beta is
#'     the parameter of the sensitivity model
#' @param beta0_fixed Optional numeric vector of values of sensitivity model
#'     intercept to profile over. If a single value, corresponds to fixing
#'     intercept at specified value. Default is NULL
#' @param weights Optional vector of patient-specific weights used for
#'     selection bias adjustment. Default is NULL
#' @param expected Whether or not to calculate the covariance matrix via the
#'     expected fisher information matrix. Default is TRUE
#' @param tol stop estimation when subsequent log-likelihood estimates are
#'     within this value
#' @param maxit Maximum number of iterations of the estimation algorithm
#' @return A "SAMBA.fit" object with nine elements: 'param', the final estimate
#' of the coeficients organized as (theta, beta),  'variance', the covariance
#' matrix of the final estimate, param.seq', the sequence of estimates at each
#' step of the EM algorithm, and 'loglik.seq', the log likelihood at each step.
#' The rest of the elements are Dstar', 'X', 'Z', and 'weights'.
#' @examples
#' library(SAMBA)
#' # These examples are generated from the vignette. See it for more details.
#'
#' # Generate IPW weights from the true model
#' expit <- function(x) exp(x) / (1 + exp(x))
#' prob.WD <- expit(-0.6 + 1 * samba.df$D + 0.5 * samba.df$W)
#' weights <- nrow(samba.df) * (1  / prob.WD) / (sum(1 / prob.WD))
#'
#' # Get initial parameter estimates
#' logit <- function(x) log(x / (1 - x))
#' fitBeta  <- glm(Dstar ~ X, binomial(), data = samba.df)
#' fitTheta <- glm(Dstar ~ Z, binomial(), data = samba.df)
#'
#' sens <- sensitivity(samba.df$Dstar, samba.df$X, mean(samba.df$D), r = 2)
#' start <- c(coef(fitTheta), logit(sens$c_marg), coef(fitBeta)[2])
#'
#' # Direct observed data likelihood maximization without fixed intercept
#' fit1 <- obsloglikEM(samba.df$Dstar, samba.df$Z, samba.df$X, start = start,
#'                  weights = weights)
#' obsloglik1 <- list(param = fit1$param, variance = diag(fit1$variance))
#'
#' # Direct observed data likelihood maximization with fixed intercept
#' fit2   <- obsloglikEM(samba.df$Dstar, samba.df$Z, samba.df$X, start = start,
#'                  beta0_fixed = logit(sens$c_marg), weights = weights)
#' # since beta0 is fixed, its variance is NA
#'
#' list(param = fit2$param, variance = diag(fit2$variance))
#' @references
#' Statistical inference for association studies using electronic health records:
#' handling both selection bias and outcome misclassification
#' Lauren J Beesley and Bhramar Mukherjee
#' medRxiv \href{https://doi.org/10.1101/2019.12.26.19015859}{2019.12.26.19015859}
#' @export
obsloglikEM <- function(Dstar, Z, X, start, beta0_fixed = NULL,
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
        if (is.data.frame(X))
            X <- as.matrix(X)
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
    theta <- start[1:(1 + ncol(Z))]
    beta  <- start[-(1:(1 + ncol(Z)))]
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

    param.seq  <- matrix(c(theta, beta), 1)
    loglik.seq <- -10 ^ 9
    while (!converged && it < maxit) {
        if (is.null(beta0_fixed)) {
            suppressWarnings({
            fit.beta <- stats::glm(Dstar ~ X, weights = p * w,
                                   family = stats::binomial())
            })
        } else {
            suppressWarnings({
            fit.beta <- stats::glm(Dstar ~ 0 + X, weights = p * w,
                                   offset = rep(beta0_fixed, length(p)),
                                   family = stats::binomial())
            })
        }
        suppressWarnings({
        fit.theta <- stats::glm(p ~ Z, family = stats::binomial(),
                                weights = weights)
        })
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

        par <- c(stats::coef(fit.theta), beta0_fixed, stats::coef(fit.beta))
        param.seq <- rbind(param.seq, par)
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

    structure(list(param = param, variance = var, param.seq = param.seq,
                   loglik.seq = loglik.seq, Dstar = Dstar, X = X, Z = Z,
                   weights = weights, beta0_fixed = beta0_fixed),
              class = "SAMBA.fit")
}
