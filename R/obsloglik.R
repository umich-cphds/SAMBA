#' Estimate parameters in the disease model using observed data log-likelihood
#' using direct maximization.
#'
#' \code{obsloglik} jointly estimates the disease model and sensitivity
#' model parameters using profile likelihood methods. Estimation involves
#' direct maximization of the observed data log-likelihood.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population of
#' interest. In this setting, we estimate parameters from the disease model
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
#' \item{Sensitivity Model}{\deqn{logit(P(D* = 1| D = 1, S = 1, X)) = beta_0 + beta_X X}}
#' }
#' @param Dstar Numeric vector containing observed disease status. Should be
#'     coded as 0/1
#' @param Z Numeric matrix of covariates in disease model
#' @param X Numeric Matrix of covariates in sensitivity model. Set to
#'     NULL to fit model with no covariates in sensitivity model
#' @param start Numeric vector of starting values for theta and beta (theta, beta).
#'     Theta is the parameter of the disease model, and beta is the parameter
#'     of the sensitivity model
#' @param beta0_fixed Optional numeric vector of values of sensitivity model
#'     intercept to profile over. If a single value, corresponds to fixing
#'     intercept at specified value. Default is NULL
#' @param weights Optional vector of patient-specific weights used for
#'     selection bias adjustment. Default is NULL
#' @param expected Whether or not to calculate the covariance matrix via the
#'     expected fisher information matrix. Default is TRUE
#' @param itnmax Maximum number of iterations to run \code{optimx}
#' @return A "SAMBA.fit" object with nine elements: 'param', the maximum
#' likelihood estimate of the coeficients,  'variance', the covariance matrix of
#' the final estimate, param.seq', the sequence of estimates at each value of
#' beta0, and 'loglik.seq', the log likelihood at each value. The rest of the
#' elements are Dstar', 'X', 'Z', and 'weights'.
#' @export
obsloglik <- function(Dstar, Z, X, start, beta0_fixed = NULL,
                          weights = NULL, expected = TRUE, itnmax = 5000)
{
    if (!is.numeric(Dstar) || !is.vector(Dstar))
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

    if (is.data.frame(Z))
        Z <- as.matrix(Z)
    if (!is.numeric(Z))
        stop("'Z' must be a numeric matrix.")
    if (is.vector(Z))
        Z <- as.matrix(Z)
    if (!is.matrix(Z))
        stop("'Z' must be a numeric matrix.")
    if (n != nrow(Z))
        stop("The number of rows of 'Z' must match the length of 'Dstar'.")

    check.weights(weights, n)
    if (is.null(weights))
        w <- rep(1, n)
    else
        w <- weights

    if (!is.logical(expected) || length(expected) > 1)
        stop("'expected' must be a length one logical.")

    if (!is.numeric(itnmax) || length(itnmax) > 1)
        stop("'itnmax' must be a single number.")

    # Define values to profile over
    if (is.null(beta0_fixed)) {
        lower  <- logit(0.2)
        upper  <- logit(0.99)
        values <- seq(lower, upper, 0.1)
    } else if (length(beta0_fixed) == 1) {
        lower  <- beta0_fixed
        upper  <- beta0_fixed
        values <- as.vector(beta0_fixed)
    } else {
        values <- as.vector(beta0_fixed)
    }

    # Loop over profile values
    param.seq  <- c()
    loglik.seq <- c()
    for (val in values) {
        #Allows us to fit intercept only model for sensitivity
        if (is.null(X)) {
            opt_long <- c(rep(NA, (1 + ncol(Z))), val)
        } else {
            opt_long <- c(rep(NA, (1 + ncol(Z))), val, rep(NA, ncol(X)))
        }

        start[2 + ncol(Z)] <- val
        opt <- optimx::optimx(start, max_obsloglik, control =
                              list(maximize = TRUE, save.failures = F, trace = 0),
                              itnmax = itnmax, method = c("BFGS", "Nelder-Mead"),
                              args = list(Z = Z, X = X, Dstar = Dstar, opt =
                                          opt_long, w = w))

        i <- which.max(opt$value)
        param.seq  <- rbind(param.seq, stats::coef(opt)[i, ])
        loglik.seq <- c(loglik.seq, opt$value[i])
    }
    param <- param.seq[which.max(loglik.seq),]

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
                   beta0_fixed = values), class = "SAMBA.fit")
}

max_obsloglik <- function(param, args)
{
    if (!is.null(args$opt))
        param[!is.na(args$opt)] <- args$opt[!is.na(args$opt)]

	Ztheta <- cbind(1, args$Z) %*% param[1:(1 + ncol(args$Z))]
	Xbeta  <- cbind(1, args$X) %*% param[-(1:(1 + ncol(args$Z)))]
    M      <- expit(Xbeta) * expit(Ztheta)

    Dstar   <- args$Dstar
    w <- args$w
    sum(w * Dstar * log(M) + w * (1 - Dstar) * log(1 - M))
}
