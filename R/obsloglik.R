#' obsloglik
#' This function can be used to jointly estimate disease model and sensitivity
#' model parameters using profile likelihood methods.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population of
#' interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
#' Notation:
#' D = binary disease status of interest
#' Dstar = observed binary disease status. Potentially a misclassified version
#' of D. We assume \code{D = 0} implies \code{Dstar = 0}.
#' S = indicator for whether patient from population of interest is included in
#'     the analytical dataset
#' Z = covariates in disease model of interest
#' W = covariates in model for patient inclusion in analytical dataset
#'     (selection model)
#' X = covariates in model for probability of observing disease given patient
#' has disease (sensitivity model)
#' Model Structure:
#' Disease Model: logit(P(D=1|X)) = theta_0 + theta_Z Z
#' Selection Model: P(S=1|W,D)
#' Sensitivity Model: logit(P(Dstar=1|D=1,X)) = beta_0 + beta_X X
#' @param Z matrix or data frame with covariates in disease model
#' @param X matrix or data frame with covariates in sensitivity model. Set to
#'     NULL to fit model with no covariates in sensitivity model.
#' @param Dstar matrix or data frame containing observed disease status
#' @param param_current vector of starting values for theta and beta (theta, beta).
#'     Theta is the parameter of the disease model, and beta is the parameter
#'     of the sensitivity model.
#' @param beta0_fixed (optional) vector of values of sensitivity model intercept
#'     to profile over. If a single value, corresponds to fixing intercept at
#'     specified value.
#' @param weights (optional) vector of subject-specific weights used for
#'     selection bias adjustment.
#' @return param vector with parameter estimates organized as (theta, beta)
#' @return param_save matrix containing estimated parameter values corresponding
#'     to different values of beta0_fixed.
#' @return loglik_save vector of log-likelihood values corresponding to
#'    different values of beta0_fixed.
#' @export

obsloglik = function(Z, X, Dstar, param_current, beta0_fixed = NULL, weights = NULL)
{
    #### Define values to profile over
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
    ### Loop over profile values
    param_save  <- c()
    loglik_save <- c()
    for (val in values) {
        #Allows us to fit intercept only model for sensitivity
        if (is.null(X)) {
            opt_long <- c(rep(NA, (1 + ncol(Z))), val)
        } else {
            opt_long <- c(rep(NA, (1 + ncol(Z))), val, rep(NA, ncol(X)))
        }

        param_current[2 + ncol(Z)] = val
        optimum <- optimx::optimx(param_current, max_obsloglik, control =
                                  list(maximize = TRUE, save.failures = F),
                                  itnmax = 5000, method = c('BFGS','Nelder-Mead'),
                                  args = list(Z = Z, X = X, Dstar = Dstar, expectedInfo =
                                  FALSE, opt = opt_long, weights = weights))

        i <- which.max(optimum$value)

        param_temp  <- as.numeric(stats::coef(optimum)[i, ])
        loglik_temp <- as.numeric(optimum$value[i])

        param_save  <- rbind(param_save, param_temp)
        loglik_save <- c(loglik_save, loglik_temp)
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
    if (is.null(weights))
        weights <- 1
	 sum(weights * Dstar * log(M) + weights * (1 - Dstar) * log(1 - M))
}
