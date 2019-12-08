#' obsloglik_var
#'
#' This function can be used to estimate the observed data
#' information matrix for joint estimation of parameters in the disease and
#' sensitivity models.
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population of
#' interest. In this setting, we estimate parameters from the disease model
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
#' @note This function assumes that the sensitivity model intercept was estimated
#' rather than fixed. If it was fixed, this will produce variances that are too big.
#' When the intercept is fixed, (1) return the information matrix, (2) remove the
#' row/column corresponding to the fixed parameter, and (3) invert the information
#' matrix to get the fixed-parameter covariance matrix
#' @param theta estimated value of theta from a call to misclass_max or
#'     misclass_maxEM
#' @param beta estimated value of beta from a call to misclass_max or
#'     misclass_maxEM
#' @param Z matrix or data frame with covariates in disease model
#' @param X matrix or data frame with covariates in sensitivity model. Set to
#'     NULL to fit model with no covariates in sensitivity model.
#' @param Dstar matrix or data frame containing observed disease status
#' @param getInfo indicator for whether the information matrix or covariance
#'     matrix should be returned. Default is FALSE.
#' @param expectedInfo indicator for whether the observed or expected
#'     information matrix should be used to estimate variance. Default is the
#'     expected information matrix (TRUE).
#'
#' @return Info estimated information matrix.
#' @return variance estimated covariance matrix.
#' @export
obsloglik_var <- function(theta, beta, X, Z, Dstar, getInfo = FALSE,
                          expectedInfo = TRUE)
{
    n <- nrow(Z)

    X1 <- cbind(1, X)
    Z1 <- cbind(1, Z)
    XBeta <- X1 %*% beta
    ZTheta <- Z1 %*% theta

    #fixed logit(1-specificity). Currently, we only support specificity = 1
    YAlpha <- matrix(-Inf, n)

    expit.zt <- expit(ZTheta)
    expit.xb <- expit(XBeta)
    expit.ya <- expit(YAlpha)

    ### Calculate Derivatives ###
    K1 <- as.vector(expit.xb * expit.zt + expit.ya * (1 - expit.zt))

    exp.xb <- exp(XBeta)
    exp.zt <- exp(ZTheta)

    dK1.dB  <- expit.xb / (1 + exp.xb) * expit.zt
    dK1.dT <- (exp.zt / ((1 + exp.zt) ^ 2)) * (expit.xb -  expit.ya)

    dK1.dBdB <- (1 - exp.xb) * exp.xb * ((1 / (1 + exp.xb)) ^ 3) * expit.zt
    dK1.dBdT <- (exp.xb / ((1 + exp.xb) ^ 2)) * (exp.zt / ((1 + exp.zt) ^ 2))
    dK1.dTdT <- (1 - exp.zt) * exp.zt * ((1 / (1 + exp.zt)) ^ 3)
    dK1.dTdT <- dK1.dTdT * (expit.xb - expit.ya)

    if (expectedInfo) {
        tmp <- as.vector(1 / (K1 * (1 - K1)))
        meat.bb <- -dK1.dB * dK1.dB * tmp
        meat.bt <- -dK1.dB * dK1.dT * tmp
        meat.tt <- -dK1.dT * dK1.dT * tmp

    } else {
        tmp <- Dstar / K1 ^ 2
        meat.bb <- tmp * (K1 * dK1.dBdB - dK1.dB * dK1.dB)
        meat.bt <- tmp * (K1 * dK1.dBdT - dK1.dB * dK1.dT)
        meat.tt <- tmp * (K1 * dK1.dTdT - dK1.dT * dK1.dT)

        tmp <- (1 - Dstar) / (1 - K1) ^ 2
        meat.bb <- meat.bb - tmp * ((1 - K1) * dK1.dBdB + dK1.dB * dK1.dB)
        meat.bt <- meat.bt - tmp * ((1 - K1) * dK1.dBdT + dK1.dB * dK1.dT)
        meat.tt <- meat.tt - tmp * ((1 - K1) * dK1.dTdT + dK1.dT * dK1.dT)
    }

    I.betabeta   <- t(apply(X1, 2, function(x1) x1 * as.vector(meat.bb))) %*% X1
    I.betatheta  <- t(apply(X1, 2, function(x1) x1 * as.vector(meat.bt))) %*% Z1
    I.thetatheta <- t(apply(Z1, 2, function(x1) x1 * as.vector(meat.tt))) %*% Z1

    info <- rbind(cbind(I.thetatheta, t(I.betatheta)),
                  cbind(I.betatheta, I.betabeta))

    if (!getInfo) {
        solve(-info)
    } else {
        -info
    }
}
