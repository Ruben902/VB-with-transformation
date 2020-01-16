rm(list=ls())
library(ggplot2)
library(INLA)
library(faraway)
library(R.matlab) 
setwd('C:/Users/rloa0001/Dropbox/VB_SkewCop_YJ/VAFC_skew/Code_Journal_Submission/Polypharm')   #mac2
DataAll = readMat('polypharm.mat')
Y = DataAll$Ytr
Y = (Y+1)/2
X = DataAll$Xtr[,1:8]
colnames(X) = c('Int','gender','race','age','M1','M2','M3','IM')
id = as.factor(matrix(t(matrix(rep(1:500,7),500,7)),7*500,1))
data = list(Y=Y,X=X,id=id)

## Priors in INLA are set on the log precision and not on the precision
## In our logistic model exp(2c) = sigma^2 = (1/tau) = exp(-logtau)
## INLA needs a prior on logtau. The prior p(c) = N(0,100) implies p(logtau) = N(0,400)
prec.prior <- list(prec = list(prior = "normal", param = c(0, 1/400)))

result  = inla(Y ~ -1+X + f(id,model="iid",hyper = prec.prior), family = c("binomial"), 
           data = data,control.family=list(link='logit'),
           control.fixed=list(mean=0, prec=(1/100), mean.intercept=0, prec.intercept=(1/100))) 

summary(result)
#Marginals for fixed effects
par(mfrow = c(3,3))
for (i in 1:length(result$marginals.fixed)) {
  tmp = inla.tmarginal(function(x) x, result$marginals.fixed[[i]]) 
  plot(tmp, type = "l", xlab = result$names.fixed[i], ylab = "Density")
}

tmp = inla.tmarginal(function(x) log(x)/(-2), result$marginals.hyperpar[[1]]) 
plot(tmp, type = "l", ylab = "Density",xlim = c(0.6,1.2),xlab = 'xi')

writeMat('INLA/polypharm_INLA.mat',PosteriorsFixed = result$marginals.fixed,PosteriorsXi = tmp)

