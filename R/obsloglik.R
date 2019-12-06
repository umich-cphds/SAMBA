#' Calculate Observed Profile Log likelihood
#'
#' \code{obsloglik} jointly estimates the disease model and sensitivity
#' model parameters using profile likelihood methods.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population of
#' interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
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
#' @param Z Numeric matrix with covariates in disease model
#' @param X Numeric Matrix with covariates in sensitivity model. Set to
#'     NULL to fit model with no covariates in sensitivity model
#' @param Dstar matrix or data frame containing observed disease status
#' @param param_current vector of starting values for theta and beta (theta, beta).
#'     Theta is the parameter of the disease model, and beta is the parameter
#'     of the sensitivity model
#' @param beta0_fixed (optional) vector of values of sensitivity model intercept
#'     to profile over. If a single value, corresponds to fixing intercept at
#'     specified value
#' @param weights (optional) vector of subject-specific weights used for
#'     selection bias adjustment
#' @param itnmax Maximum number of iterations to run \code{optimx}
#' @return param vector with parameter estimates organized as (theta, beta)
#' @return param_save matrix containing estimated parameter values corresponding
#'     to different values of beta0_fixed.
#' @return loglik_save vector of log-likelihood values corresponding to
#'    different values of beta0_fixed.
#' @export
obsloglik <- function(Dstar, Z, X, param_current, weights = NULL,
                      beta0_fixed = NULL, itnmax = 5000)
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
        weights <- rep(1, n)

    if (!is.numeric(itnmax) || length(itnmax) > 1)
        stop("'itnmax' must be a single number.")

    # Define values to profile over
    if (is.null(beta0_fixed)) {
        lower  <- logit(0.2)
        upper  <- logit(0.99)
        values <- as.vector(seq(lower, upper, 0.1))
    } else if (length(beta0_fixed) == 1) {
        lower  <- beta0_fixed
        upper  <- beta0_fixed
        values <- as.vector(beta0_fixed)
    } else {
        values <- as.vector(beta0_fixed)
    }

    # Loop over profile values
    param_save  <- c()
    loglik_save <- c()
    for (val in values) {
        #Allows us to fit intercept only model for sensitivity
        if (is.null(X)) {
            opt_long <- c(rep(NA, (1 + ncol(Z))), val)
        } else {
            opt_long <- c(rep(NA, (1 + ncol(Z))), val, rep(NA, ncol(X)))
        }

        param_current[2 + ncol(Z)] <- val
        opt <- optimx::optimx(param_current, max_obsloglik, control =
                              list(maximize = TRUE, save.failures = F), itnmax =
                              itnmax, method = c("BFGS", "Nelder-Mead"),
                              args = list(Z = Z, X = X, Dstar = Dstar,
                                          opt = opt_long, weights = weights))

        i <- which.max(opt$value)
        param_save  <- rbind(param_save, stats::coef(opt)[i, ])
        loglik_save <- c(loglik_save, opt$value[i])
    }
    param_new = param_save[which.max(loglik_save),]
    list(param = param_new, param_save = param_save, loglik_save = loglik_save)
}

max_obsloglik <- function(param, args)
{
    if (!is.null(args$opt))
        param[!is.na(args$opt)] <- args$opt[!is.na(args$opt)]

	Ztheta <- cbind(1, args$Z) %*% param[1:(1 + ncol(args$Z))]
	Xbeta  <- cbind(1, args$X) %*% param[-(1:(1 + ncol(args$Z)))]
    M      <- expit(Xbeta) * expit(Ztheta)

    Dstar   <- args$Dstar
    weights <- args$weights
    sum(weights * Dstar * log(M) + weights * (1 - Dstar) * log(1 - M))
}
