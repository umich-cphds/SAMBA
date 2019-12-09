#' obsloglik_weightedvar
#'
#' This function can be used to estimate the observed data information matrix
#' for joint estimation of parameters in the disease and sensitivity models.
#'
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population
#' of interest. In this setting, we estimate parameters from the disease model
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
#'
#' @note This function assumes that the sensitivity model intercept was
#' estimated rather than fixed. If it was fixed, this will produce variances
#' that are too big. When the intercept is fixed, (1) return the information
#' matrix, (2) remove the row/column corresponding to the fixed parameter,
#' and (3) invert the information matrix to get the fixed-parameter covariance
#' matrix
#'
#' @param Dstar Numeric vector containing observed disease status. Should be
#'     coded as 0/1
#' @param Z Numeric matrix of covariates in disease model
#' @param X Numeric matrix with covariates in sensitivity model. Set to NULL
#'     to fit model with no covariates in sensitivity model. 'X' should not
#'     contain an intercept
#' @param theta estimated value of theta from a call to misclass_max or misclass_maxEM
#' @param beta estimated value of beta from a call to misclass_max or misclass_maxEM
#' @param weights Optional vector of subject-specific weights used for
#'     selection bias adjustment. Default is NULL
#' @param getInfo Logical variable for whether the information matrix or
#'     covariance matrix should be returned. Default is FALSE.
#' @param expectedInfo Logical variable for whether the observed or expected
#'     information matrix should be used to estimate variance. Default is the
#'     expected information matrix (TRUE).
#'
#' @return Info estimated information matrix.
#' @return variance estimated covariance matrix.
#' @export
obsloglik_weightedvar <- function(Dstar, Z, X, theta, beta, weights = NULL,
                                      getInfo = FALSE, expectedInfo = TRUE)
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

    if (!is.logical(getInfo))
        stop("'getInfo' must take a logical value.")

    if (!is.logical(expectedInfo))
        stop("'expectedInfo' must take a logical value.")

    X1 <- cbind(1, X)
    Z1 <- cbind(1, Z)

    XBeta  <- X1 %*% beta
    ZTheta <- Z1 %*% theta

    if (!is.numeric(beta) || !is.vector(beta))
        stop("'beta' must be a numeric vector.")
    if (length(beta) != ncol(X1))
        stop("'beta' must have length 'ncol(X) + 1'.")

    if (!is.numeric(theta) || !is.vector(theta))
        stop("'theta' must be a numeric vector.")
    if (length(theta) != ncol(Z1))
        stop("'theta' must have length 'ncol(Z) + 1'.")

    # fixed logit(1-specificity). Currently, we only support specificity = 1
    YAlpha <- matrix(-Inf, n)

    expit.xb <- expit(XBeta)
    expit.zt <- expit(ZTheta)
    expit.ya <- expit(YAlpha)

    exp.xb <- exp(XBeta)
    exp.zt <- exp(ZTheta)

    K1 <- as.vector(expit.xb * expit.zt + expit.ya * (1 - expit.zt))

    dK1.dB <- expit.xb * (1 / (1 + exp.xb)) * expit.zt
    dK1.dT <- exp.zt / (1 + exp.zt) ^ 2 * (expit.xb - expit.ya)

    dK1.dBdB <- (1 - exp.xb) * exp.xb / (1 + exp.xb) ^ 3 * expit.zt
    dK1.dBdT <- exp.xb / (1 + exp.xb) ^ 2 * exp.zt / (1 + exp.zt) ^ 2
    dK1.dTdT <- (1 - exp.zt) * exp.zt / (1 + exp.zt) ^ 3 * (expit.xb - expit.ya)

    # Calculate information matrix
    if (expectedInfo) {
        tmp <- 1 / (K1 * (1 - K1))
        meat.bb <- -dK1.dB * dK1.dB * tmp
        meat.bt <- -dK1.dB * dK1.dT * tmp
        meat.tt <- -dK1.dT * dK1.dT * tmp
    } else {
        tmp <- 1 / K1 ^ 2
        meat.bb <- Dstar * (K1 * dK1.dBdB - dK1.dB * dK1.dB) * tmp
        meat.bt <- Dstar * (K1 * dK1.dBdT - dK1.dB * dK1.dT) * tmp
        meat.tt <- Dstar * (K1 * dK1.dTdT - dK1.dT * dK1.dT) * tmp

        tmp <- 1 / (1 - K1) ^ 2
        meat.bb <- meat.bb - (1 - Dstar) * (1 - K1) * dK1.dBdB +
                       dK1.dB * dK1.dB * tmp
        meat.bt <- meat.bt - (1 - Dstar) * (1 - K1) * dK1.dBdT +
                       dK1.dB * dK1.dT * tmp
        meat.tt <- meat.tt - (1 - Dstar) * (1 - K1) * dK1.dTdT +
                       dK1.dT * dK1.dT * tmp
    }

    # BREAD
    I.betabeta <- t(apply(X1, 2, function(x) x * weights * meat.bb)) %*% X1
    I.betatheta <- t(apply(X1, 2, function(x) x * weights * meat.bt)) %*% Z1
    I.thetatheta <- t(apply(Z1, 2, function(x) x * weights * meat.tt)) %*% Z1
    info <- rbind(cbind(I.thetatheta, t(I.betatheta)),
                  cbind(I.betatheta, I.betabeta))

    # MEAT
    tmp <- (Dstar - K1) / (K1 * (1 - K1))
    U.beta <- apply(X1, 2, function(x) tmp * dK1.dB * x)
    U.theta <- apply(Z1, 2, function(z) tmp * dK1.dT * z)
    U <- cbind(U.theta, U.beta)

    meat <- t(U) %*% diag(weights ^ 2) %*% U

    if (!getInfo) {
        solve(-info) %*% meat %*% solve(-info)
    } else {
        list(inv_bread = -info, meat = meat)
    }
}
