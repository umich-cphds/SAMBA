
#' sensitivity
#' @description This function can be used to estimate the marginal sensitivity and sensitivity as a function of covariates X for a misclassified binary outcome.
#'
#' @param X matrix or data frame with covariates in sensitivity model. Set to NULL to fit model with no covariates in sensitivity model.
#' @param Dstar matrix or data frame containing observed disease status
#' @param r (optional) marginal sampling ratio, P( sampled | disease )/P( sampled | no disease ). Only one of r and weights can be specified.
#' @param weights (optional) vector of subject-specific weights used for selection bias adjustment. Only one of r and weights can be specified.
#' @param prev disease prevalence in the population or subject-specific P( disease | X ) in population
#'
#' @return c_marg marginal sensitivity estimate, P( observe disease | have disease )
#' @return c_X sensitivity as a function of X, P( observe disease | have disease, X )
#'
#' @details We are interested in modeling the relationship between binary disease status and covariates Z using a logistic regression model. However, D may be misclassified, and our observed data may not well-represent the population of interest. In this setting, we estimate parameters from the disease model using the following modeling framework.
#' @details Notation:
#' @details D = binary disease status of interest
#' @details Dstar = observed binary disease status. Potentially a misclassified version of D. We assume D=0 implies Dstar=0.
#' @details S = indicator for whether patient from population of interest is included in the analytical dataset
#' @details Z = covariates in disease model of interest
#' @details W = covariates in model for patient inclusion in analytical dataset (selection model)
#' @details X = covariates in model for probability of observing disease given patient has disease (sensitivity model)
#' @details Model Structure:
#' @details Disease Model: logit(P(D=1|X)) = theta_0 + theta_Z Z
#' @details Selection Model: P(S=1|W,D)
#' @details Sensitivity Model: logit(P(Dstar=1|D=1,X)) = beta_0 + beta_X X
#'
#' @export
sensitivity <- function(X, Dstar, r = NULL, weights = NULL, prev) {
  myLinkprev = modLinkprev()
  if(length(prev)>1){
    print('Using average prevalence to calculate marginal c')
  }
  if(is.null(r) & is.null(weights)){
    p_star = mean(Dstar)
    c_marg = p_star/mean(prev)
    c_marg = ifelse(c_marg > 1, 1, c_marg)
    fitBeta = glm(Dstar~as.matrix(X), family = binomial())
    starting  = c(logit(c_marg), coef(fitBeta)[2:length(coef(fitBeta))])
    prevr <<- prev
    fitSENS = try(glm(Dstar ~ as.matrix(X),family=binomial(link=myLinkprev), start = starting),silent = TRUE)
    if(class(fitSENS)[1] != "try-error"){
      c1 = expit(cbind(1,X) %*% as.matrix(coef(fitSENS)))
    }else{
      c1 = predict(fitBeta, type = 'response')/prev
    }
  }else if(is.null(r) & !is.null(weights)){
    p_star = sum(Dstar*weights)/sum(weights)
    c_marg = p_star/mean(prev)
    c_marg = ifelse(c_marg > 1, 1, c_marg)
    fitBeta = glm(Dstar~as.matrix(X), family = binomial(), weights = weights)
    starting  = c(logit(c_marg),  coef(fitBeta)[2:length(coef(fitBeta))])
    prevr <<- prev
    fitSENS = try(glm(Dstar ~ as.matrix(X),family=binomial(link=myLinkprev), start = starting, weights = weights),silent = TRUE)
    if(class(fitSENS)[1] != "try-error"){
      c1 = expit(cbind(1,X) %*% as.matrix(coef(fitSENS)))
    }else{
      c1 = predict(fitBeta, type = 'response')/prev
    }
    c1 = ifelse(c1 > 1, rep(1, length(Dstar)), c1)
  }else if(!is.null(r) & is.null(weights)){
    p_star = mean(Dstar)
    prevr = r*prev/(r*prev + 1-prev)
    c_marg = p_star/mean(prevr)
    c_marg = ifelse(c_marg > 1, 1, c_marg)
    fitBeta = glm(Dstar~as.matrix(X), family = binomial())
    starting  = c(logit(c_marg),  coef(fitBeta)[2:length(coef(fitBeta))])
    prevr <<- r*prev/(r*prev + 1-prev)
    fitSENS = try(glm(Dstar ~ as.matrix(X),family=binomial(link=myLinkprev), start = starting),silent = TRUE)
    if(class(fitSENS)[1] != "try-error"){
      c1 = expit(cbind(1,as.matrix(X)) %*% as.matrix(coef(fitSENS)))
    }else{
      c1 = predict(fitBeta, type = 'response')/prevr
    }
  }else{
    stop('Error: only one of r and weights can be specified')
  }
  c_marg = ifelse(c_marg < 1, c_marg, NA)
  c1 = ifelse(c1 > 1, rep(1, length(Dstar)), c1)
  return(list(c_marg = c_marg, c_X = c1))
}
