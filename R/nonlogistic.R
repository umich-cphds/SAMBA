#' Estimate parameters in the disease model given sensitivity as a function of
#' covariates.
#'
#' non-logistic link function for D* given Z and sensitivity. This function
#' assumes that sensitivity as a function of X is known or has been estimated
#'
#' We are interested in modeling the relationship between binary disease status
#' and covariates Z using a logistic regression model. However, D may be
#' misclassified, and our observed data may not well-represent the population
#' of interest. In this setting, we estimate parameters from the disease model
#' using the following modeling framework.
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
#' @param Z numeric matrix of covariates in disease model
#' @param weights Optional numeric vector of patient-specific weights used for
#'     selection bias adjustment. Default is NULL
#' @param c_X sensitivity as a function of X, P(D* = 1| D = 1, S = 1, X)
#' @return a list with two elements: (1) 'param', a vector with parameter
#'     estimates for disease model (logOR of Z), and (2) 'variance', a vector of
#'     variance estimates for disease model parameters. Results do not include
#'     intercept.
#' @examples
#' library(SAMBA)
#' # These examples are generated from the vignette. See it for more details.
#'
#' # Generate IPW weights from the true model
#' expit <- function(x) exp(x) / (1 + exp(x))
#' prob.WD <- expit(-0.6 + 1 * samba.df$D + 0.5 * samba.df$W)
#' weights <- nrow(samba.df) * (1  / prob.WD) / (sum(1 / prob.WD))
#'
#' # Estimate sensitivity by using inverse probability of selection weights
#' # and P(D=1)
#' sens <- sensitivity(samba.df$Dstar, samba.df$X, prev = mean(samba.df$D),
#'                     weights = weights)
#'
#' nonlog1 <- nonlogistic(samba.df$Dstar, samba.df$Z, c_X = sens$c_X,
#'                        weights = weights)
#' @references
#' Statistical inference for association studies using electronic health records:
#' handling both selection bias and outcome misclassification
#' Lauren J Beesley, Bhramar Mukherjee
#' medRxiv \href{https://doi.org/10.1101/2019.12.26.19015859}{2019.12.26.19015859}
#' @export
nonlogistic <- function(Dstar, Z, c_X, weights = NULL)
{
    if (!is.numeric(Dstar) || !is.vector(Dstar))
        stop("'Dstar' must be a numeric vector.")
    if (length(setdiff(0:1, unique(Dstar))) != 0)
        stop("'Dstar' must be coded 0/1.")

    n <- length(Dstar)
    if (is.data.frame(Z))
        Z <- as.matrix(Z)
    if (!is.numeric(Z))
        stop("'Z' must be a numeric matrix.")
    if (is.vector(Z))
        Z <- as.matrix(Z)
    if (!is.matrix(Z))
        stop("'Z' must be a numeric matrix.")

    if (nrow(Z) != n)
        stop("The number of rows of 'Z' must match the length of 'Dstar'.")

    if (!is.numeric(c_X) || !is.vector(c_X))
        stop("'c_X' must be a numeric vector.")
    if (length(c_X) != length(Dstar))
        stop("'c_X' must have the same length as 'Dstar'.")

    check.weights(weights, n)

    fitTheta <- stats::glm(Dstar ~ Z, family = stats::binomial())
    start <- stats::coef(fitTheta)
    if (is.null(weights)) {
        fit <- stats::glm(Dstar ~ Z, start = start, family =
                          stats::binomial(modLink(c_X)))

        var <- diag(summary(fit)$cov.scaled)
    } else {
        # Estimate using survey package if weighted
        colnames(Z) <- paste0("Z", 1:ncol(Z))
        design <- survey::svydesign(ids = 1:length(Dstar), strata = NULL,
                                    data = data.frame(Z, Dstar = Dstar),
                                    weights = weights)

        form <- paste("Dstar ~", paste(colnames(Z), collapse = "+"))
        fit  <- survey::svyglm(stats::formula(form), design = design,
                               family = stats::binomial(modLink(c_X)),
                               start = start)

        var <- diag(stats::vcov(fit))
    }

    param <- stats::coef(fit)
    list(param = param, variance = var)
}
