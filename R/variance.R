obsloglik_var <- function(Dstar, Z, X, theta, beta, beta0_fixed,
                              expected = TRUE)
{
    fixed <- !is.null(beta0_fixed) & length(beta0_fixed) < 2

    n <- length(Dstar)

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

    if (expected) {
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

    if (fixed) {
        k <- ncol(Z) + 2
        var <- -solve(info[-k , -k])
        var <- cbind(var[, 1:(k - 1)], NA, var[, k:ncol(var)])
        var <- rbind(var[1:(k - 1),], NA, var[k:nrow(var),])
    } else {
        var <- -solve(info)
    }

    var
}

obsloglik_var_weighted <- function(Dstar, Z, X, theta, beta, beta0_fixed,
                                       weights, expected = TRUE)
{

    if (!is.logical(expected) || length(expected) > 1)
        stop("'expected' must be a length one logical.")

    fixed <- !is.null(beta0_fixed) && length(beta0_fixed) == 1

    X1 <- cbind(1, X)
    Z1 <- cbind(1, Z)

    XBeta <- X1 %*% beta
    ZTheta <- Z1 %*% theta

    n <- length(Dstar)
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
    if (expected) {
        tmp <- 1 / (K1 * (1 - K1))
        bread.bb <- -dK1.dB * dK1.dB * tmp
        bread.bt <- -dK1.dB * dK1.dT * tmp
        bread.tt <- -dK1.dT * dK1.dT * tmp
    } else {
        tmp <- 1 / K1 ^ 2
        bread.bb <- Dstar * (K1 * dK1.dBdB - dK1.dB * dK1.dB) * tmp
        bread.bt <- Dstar * (K1 * dK1.dBdT - dK1.dB * dK1.dT) * tmp
        bread.tt <- Dstar * (K1 * dK1.dTdT - dK1.dT * dK1.dT) * tmp

        tmp <- 1 / (1 - K1) ^ 2
        bread.bb <- bread.bb - (1 - Dstar) * (1 - K1) * dK1.dBdB +
                       dK1.dB * dK1.dB * tmp
        bread.bt <- bread.bt - (1 - Dstar) * (1 - K1) * dK1.dBdT +
                       dK1.dB * dK1.dT * tmp
        bread.tt <- bread.tt - (1 - Dstar) * (1 - K1) * dK1.dTdT +
                       dK1.dT * dK1.dT * tmp
    }

    # BREAD
    I.betabeta <- t(apply(X1, 2, function(x) x * weights * bread.bb)) %*% X1
    I.betatheta <- t(apply(X1, 2, function(x) x * weights * bread.bt)) %*% Z1
    I.thetatheta <- t(apply(Z1, 2, function(x) x * weights * bread.tt)) %*% Z1
    info <- rbind(cbind(I.thetatheta, t(I.betatheta)),
                  cbind(I.betatheta, I.betabeta))

    # MEAT
    tmp <- (Dstar - K1) / (K1 * (1 - K1))
    U.beta <- apply(X1, 2, function(x) tmp * dK1.dB * x)
    U.theta <- apply(Z1, 2, function(z) tmp * dK1.dT * z)
    U <- cbind(U.theta, U.beta)

    U <- apply(U, 2, function(u) u * weights)

    meat <- t(U) %*% U

    var <- solve(-info) %*% meat %*% solve(-info)

    if (fixed) {
        k <- ncol(Z) + 2
        bread <- -solve(info[-k, -k])
        var <- bread %*% meat[-k, -k] %*% bread
        var <- cbind(var[, 1:(k - 1)], NA, var[, k:ncol(var)])
        var <- rbind(var[1:(k - 1),], NA, var[k:nrow(var),])
    } else {
        var <- -solve(info) %*% meat %*% -solve(info)
    }
    var
}
