## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, fig.cap="Figure 1: Model Structure", out.width = '80%'----
knitr::include_graphics("images/ModelDiagram.png")

## ---- echo = TRUE, eval = TRUE,  fig.width = 7, fig.height= 4------------
library(SAMBA)
options(warn = -1)
mSample = function(Y){return(sample(x=c(0,1), size = 1, prob = c(1-Y,Y)))}
expit = function(x){return(exp(x)/(1+exp(x)))}
logit = function(x){log(x/(1-x))}
Nobs = 5000

### Generate Predictors and Follow-up Information
set.seed('1234')
COV = MASS::mvrnorm(n=Nobs, mu = c(0,0,0), Sigma = rbind(c(1,0,0.4), c(0,1,0), c(0.4,0,1)))
DAT = data.frame(Z = COV[,1], X = COV[,2], W = COV[,3])

### Generate random uniforms
set.seed('5678'); U1 = runif(n = Nobs, min = 0, max = 1);
set.seed('4321'); U2 = runif(n = Nobs, min = 0, max = 1);
set.seed('8765'); U3 = runif(n = Nobs, min = 0, max = 1);

### Generate Disease Status
DISEASE = expit(-2+0.5*DAT$Z)
DAT$D = ifelse(DISEASE > U1, 1, 0) #apply(matrix(probD),1, mSample)

### Relate W and D
DAT$W = DAT$W + 1*DAT$D

### Generate Misclassification
SENS = expit(-0.4 + 1*DAT$X)
SENS[DAT$D == 0] = 0
DAT$Dstar = ifelse(SENS > U2, 1, 0)# apply(matrix(prob1),1, mSample)

### Generate Sampling Status
SELECT = expit(-0.6 +1*DAT$D + 0.5*DAT$W )
DAT$S = ifelse(SELECT > U3, 1, 0) #apply(matrix(prob1),1, mSample)

### Observed Data
DAT_SUB = DAT[DAT$S==1,]

### True marginal sampling ratio
prob1 = expit(-0.6 +1*1 + 0.5*DAT$W )
prob0 = expit(-0.6 +1*0 + 0.5*DAT$W )
r_marg_true = mean(prob1[DAT$D==1])/mean(prob0[DAT$D==0]) 

### True inverse probability of sampling weights
prob_WD =expit(-0.6 +1*DAT_SUB$D + 0.5*DAT_SUB$W )
weights = length(DAT_SUB[,1])*(1/prob_WD)/(sum(1/prob_WD))

### True associations with D in population
trueX = glm(D~X, family = binomial(), data = DAT)
trueZ = glm(D~Z, family = binomial(), data = DAT)

### Initial Parameter Values
fitBeta = glm(Dstar~X, family = binomial(),data = DAT_SUB)
fitTheta = glm(Dstar~Z, family = binomial(), data = DAT_SUB)

## ---- echo = TRUE, eval = TRUE,  fig.width = 7, fig.height= 4,  results='hide'----
### Using marginal sampling ratio r and P(D=1)
sens1 = sensitivity(Z = DAT_SUB$Z, X = DAT_SUB$X, Dstar = DAT_SUB$Dstar, 
                    r = r_marg_true, prev = mean(DAT$D))
### Using inverse probability of selection weights and P(D=1)
sens2 = sensitivity(Z = DAT_SUB$Z, X = DAT_SUB$X, Dstar = DAT_SUB$Dstar, 
                    weights = weights, prev = mean(DAT$D))
### Using marginal sampling ratio r and P(D=1|X)
sens3 = sensitivity(Z = DAT_SUB$Z, X = DAT_SUB$X, Dstar = DAT_SUB$Dstar, 
                    r = r_marg_true, prev = predict(trueX, newdata = DAT_SUB, type = 'response'))
### Using inverse probability of selection weights and P(D=1|X)
sens4 = sensitivity(Z = DAT_SUB$Z, X = DAT_SUB$X, Dstar = DAT_SUB$Dstar, 
                    weights = weights, prev = predict(trueX, newdata = DAT_SUB, type = 'response'))

## ---- echo = TRUE, eval = TRUE,  fig.width = 7, fig.height= 4, message=FALSE, results='hide'----
### Approximation of D*|Z
approx1 = approxdist(Z = DAT_SUB$Z, Dstar = DAT_SUB$Dstar, weights = weights, c_marg = sens1$c_marg)

### Non-logistic link function method
nonlog1 = nonlogistic(Z = DAT_SUB$Z, X = DAT_SUB$X, Dstar = DAT_SUB$Dstar, weights = weights, c_X = sens3$c_X)

### Direct observed data likelihood maximization with fixed intercept
starting = as.numeric(c(coef(fitTheta), logit(sens1$c_marg),coef(fitBeta)[2]))
fit = obsloglik(Z=DAT_SUB$Z, X=DAT_SUB$X, Dstar  = DAT_SUB$Dstar, 
                   param_current = starting,
                   beta0_fixed = logit(sens1$c_marg), weights = weights)
Info = obsloglik_var(theta = fit$param[1:2], beta = fit$param[3:4],
                    X=DAT_SUB$X, Z = DAT_SUB$Z, Dstar = DAT_SUB$Dstar, getInfo = TRUE, weights = weights)
obsloglik1 = list(param = fit$param, variance = diag(solve(-Info[c(1,2,4),c(1,2,4)])))

## ---- echo = FALSE, eval = TRUE,  fig.width = 5, fig.height= 5-----------
plot(sort(sens3$c_X), xlab = 'Patients', ylab = 'Sensitivity', main = 'Figure 2: Sensitivity Estimates', type = 'l', col = 'red', lwd = 2)
lines(sort(expit(obsloglik1$param[3] + obsloglik1$param[4]*DAT_SUB$X)), col = 'blue', lwd = 2)
abline(h=sens1$c_marg, col = 'purple', lwd = 2)
lines(sort(expit(-0.4 + 1*DAT_SUB$X)), col = 'black', lwd = 2)
legend(x='topleft', fill = c('purple', 'red','blue', 'black'), legend = c('Estimated marginal sensitivity', 'Estimated using non-logistic link method', 'Estimated using observed data log-likelihood', 'Truth'), cex = 0.8)

## ---- echo = FALSE, eval = TRUE,  fig.width = 5, fig.height= 5,  results='hide'----
rvals = c(1,1.5,2,2.5,5,10)
COL = c('red', 'orange', 'yellow', 'green', 'blue', 'purple')
true_prevs = predict(trueX, newdata = DAT_SUB, type = 'response')
plot(sort(expit(-0.4 + 1*DAT_SUB$X)), xlab = 'Patients', ylab = 'Sensitivity', main = 'Figure 3: Estimated sensitivity across \n marginal sampling ratios', type = 'l', col = 'black', lwd = 2, ylim = c(0,1))
for(i in 1:length(rvals)){
  TEMP = sensitivity(Z = DAT_SUB$Z, X = DAT_SUB$X, Dstar = DAT_SUB$Dstar,  r = rvals[i], prev = true_prevs)
  lines(sort(TEMP$c_X), col = COL[i])}
legend(x='topleft', legend = c(rvals, 'Truth'), title = 'Sampling Ratio', fill = c(COL, 'black'), cex = 0.8)

## ---- echo = FALSE, eval = TRUE,  fig.width = 7, fig.height= 4, message=FALSE----
library(ggplot2)
library(scales)

## ---- echo = FALSE, eval = TRUE,  fig.width = 7, fig.height= 4-----------
METHODS = c('True',  'Uncorrected', 'Approx D*|Z + IPW','Non-logistic Link + IPW','Obs. log-lik + IPW')
PARAM = c( coef(trueZ)[2], coef(fitTheta)[2],approx1$param,  nonlog1$param[2], obsloglik1$param[2] )
VARIANCE = c(diag(summary(trueZ)$cov.scaled)[2],diag(summary(fitTheta)$cov.scaled)[2],
             approx1$variance,nonlog1$variance[2],obsloglik1$variance[2])
pd = position_dodge(width=0.6)
a <- ggplot(data = data.frame(METHODS = METHODS, PARAM = PARAM, VARIANCE = VARIANCE),
       aes(xmin= METHODS, xmax = METHODS, ymin = PARAM - 1.96*sqrt(VARIANCE), ymax =  PARAM + 1.96*sqrt(VARIANCE),
           col = METHODS,x = METHODS, y = PARAM)) +
  geom_point(position = position_dodge(.7), size = 2) +
  geom_linerange(position = position_dodge(.7), size = 1.2) +
  xlab('') + ylab('logOR')+ggtitle('Figure 4: Estimated Log-Odds Ratio Across Methods')+
  scale_x_discrete(limits=METHODS)+
  geom_hline(yintercept = PARAM[1], linetype = 1, color = 'black')+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme(legend.position="top",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.text=element_text(size=8), legend.title = element_blank(),
        axis.text.x=element_text(angle=20,hjust=1,vjust=1), text = element_text(size=12))
print(a)

