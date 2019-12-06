expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1 - x))

modLinkprev <- function(prevr)
{
    linkfun <- function(y) log(y / (prevr-y))
    linkinv <- function(eta) prevr * exp(eta) / (1 + exp(eta))

    # derivative of invlink wrt eta
    mu.eta   <- function(eta) prevr * exp(eta) / (1 + exp(eta))^2
    valideta <- function(eta) TRUE

    link <- "log(y/(p(D=1|X,S=1)-y))"
    structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                   valideta = valideta, name = link),
              class = "link-glm")
}

modLink <- function(cX)
{
    linkfun <- function(y) log(y / (cX - y))
    linkinv <- function(eta) cX * exp(eta) / (1 + exp(eta))

    # derivative of invlink wrt eta
    mu.eta   <- function(eta) cX * exp(eta) / (1 + exp(eta)) ^ 2
    valideta <- function(eta) TRUE

    link <- "log(y/(c_X-y))"
    structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                   valideta = valideta, name = link),
              class = "link-glm")
}

check.weights <- function(weights, n)
{
    if (!is.null(weights)) {
        if (length(weights) != n)
            stop("The length of 'weights' must match the length of 'Dstar'.")
        if (!is.numeric(weights) || !is.vector(weights))
            stop("'weights' must be a numeric vector.")
        if (any(weights < 0))
            stop("'weights' must be nonnegative.")
    }
}
