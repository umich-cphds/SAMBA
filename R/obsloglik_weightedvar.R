#' obsloglik_weightedvar
#' @description This function can be used to estimate the observed data information matrix for joint estimation of parameters in the disease and sensitivity models.
#'
#' @param theta estimated value of theta from a call to misclass_max or misclass_maxEM
#' @param beta estimated value of beta from a call to misclass_max or misclass_maxEM
#' @param Z matrix or data frame with covariates in disease model
#' @param X matrix or data frame with covariates in sensitivity model. Set to NULL to fit model with no covariates in sensitivity model.
#' @param Dstar matrix or data frame containing observed disease status
#' @param getInfo indicator for whether the information matrix or covariance matrix should be returned. Default is FALSE.
#' @param expectedInfo indicator for whether the observed or expected information matrix should be used to estimate variance. Default is the expected information matrix (TRUE).
#' @param weights (optional) vector of subject-specific weights used for selection bias adjustment.
#'
#' @return Info estimated information matrix.
#' @return variance estimated covariance matrix.
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
#' @details This function assumes that the sensitivity model intercept was estimated rather than fixed. If it was fixed, this will produce variances that are too big. When the intercept is fixed, (1) return the information matrix, (2) remove the row/column corresponding to the fixed parameter, and (3) invert the information matrix to get the fixed-parameter covariance matrix
#'
#' @export

obsloglik_weightedvar = function(theta, beta, X,Z,Dstar, getInfo = FALSE, expectedInfo = TRUE, weights = NULL){
  Nobs = length(as.matrix(Z)[,1])
  XBeta = as.matrix(cbind(rep(1,Nobs),X)) %*% as.matrix(beta)

  alpha_fixed = -Inf #fixed logit(1-specificity). Currently, we only support specificity = 1
  YAlpha = cbind(as.matrix(rep(alpha_fixed, Nobs)))
  ZTheta = as.matrix(cbind(rep(1,Nobs),Z)) %*% as.matrix(theta)

  #############################
  ### Calculate Derivatives ###
  #############################
  K1 = expit(XBeta)*expit(ZTheta)+expit(YAlpha)*(1-expit(ZTheta))

  dK1_beta = as.vector(expit(XBeta)*(1/(1+exp(XBeta)))*expit(ZTheta))#check
  dK1_theta = as.vector(     (exp(ZTheta)/((1+exp(ZTheta))^2)) * ( expit(XBeta) -  expit(YAlpha)) ) #check

  dK1_betabeta = as.vector((1-exp(XBeta))*exp(XBeta)*((1/(1+exp(XBeta)))^3)  *expit(ZTheta))
  dK1_betatheta = as.vector(  (exp(XBeta)/((1+exp(XBeta))^2))   *     (exp(ZTheta)/((1+exp(ZTheta))^2)))
  dK1_thetatheta = as.vector(  (1-exp(ZTheta))*exp(ZTheta)*((1/(1+exp(ZTheta)))^3)  *  (expit(XBeta) - expit(YAlpha)))


  #############################################
  ### Calculate Expected Information Matrix ###
  #############################################
  if(expectedInfo){
    meat_betabeta = -as.vector(dK1_beta * dK1_beta) * as.vector(1 / (K1 * (1 - K1)))
    meat_betatheta =-as.vector(dK1_beta* dK1_theta)*as.vector(1/(K1*(1-K1)))
    meat_thetatheta =-as.vector(dK1_theta* dK1_theta)*as.vector(1/(K1*(1-K1)))
    #exp(XBeta+ZTheta)/((1+exp(ZTheta))*(1+exp(XBeta))*(1+exp(XBeta)+exp(ZTheta)))
  }else{
    meat_betabeta =Dstar * (( K1*dK1_betabeta - dK1_beta*dK1_beta   )/(K1^2))
    meat_betatheta =Dstar * (( K1*dK1_betatheta - dK1_beta*dK1_theta   )/(K1^2))
    meat_thetatheta =Dstar * (( K1*dK1_thetatheta - dK1_theta*dK1_theta   )/(K1^2))
    meat_betabeta =meat_betabeta - (1-Dstar) * (( (1-K1)*dK1_betabeta + dK1_beta*dK1_beta   )/((1-K1)^2))
    meat_betatheta =meat_betatheta - (1-Dstar) * (( (1-K1)*dK1_betatheta + dK1_beta*dK1_theta   )/((1-K1)^2))
    meat_thetatheta =meat_thetatheta - (1-Dstar) * (( (1-K1)*dK1_thetatheta + dK1_theta*dK1_theta   )/((1-K1)^2))
  }

  ### BREAD
   I_betabeta = t( sweep(as.matrix(cbind(rep(1,Nobs),X)), MARGIN=1, weights*meat_betabeta, `*`) ) %*% as.matrix(cbind(rep(1,Nobs),X))
   I_betatheta = t(sweep(as.matrix(cbind(rep(1,Nobs),X)), MARGIN=1, weights*meat_betatheta, `*`))  %*% as.matrix(cbind(rep(1,Nobs),Z))
   I_thetatheta = t( sweep(as.matrix(cbind(rep(1,Nobs),Z)), MARGIN=1, weights*meat_thetatheta, `*`)) %*% as.matrix(cbind(rep(1,Nobs),Z))
   Info = rbind(cbind(I_thetatheta, t(I_betatheta)),
                cbind(I_betatheta, I_betabeta))

  ### MEAT
  U_beta = sweep(as.matrix(cbind(rep(1,Nobs),X)), MARGIN = 1, ((Dstar - K1)/(K1*(1-K1)))*dK1_beta, `*`)
  U_theta =sweep(as.matrix(cbind(rep(1,Nobs),Z)), MARGIN = 1, ((Dstar - K1)/(K1*(1-K1)))*dK1_theta, `*`)
  U = cbind(U_theta, U_beta)

  prodapply = function(i){
    return( as.matrix(U[i,]*weights[i]) %*% t(as.matrix(U[i,]*weights[i])))
  }
  products = lapply(as.matrix(1:Nobs),prodapply)
  prodsum = apply(simplify2array(products),1:2,sum)

   if(!getInfo){
    variance = solve(-Info) %*% prodsum %*% solve(-Info)
    return(variance)
  }else{
    return(list(inv_bread = as.matrix(-Info),meat =prodsum  ))
  }
}
