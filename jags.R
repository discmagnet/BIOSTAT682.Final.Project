# JAGS Implementation of the AFT Model
library(R2jags)
library(survival)
library(dplyr)
pbc<-read.csv("pbc.csv")

pbc$drug <- 1*(pbc$treatment==1)
pbc$female <- 1*(pbc$sex == 1)
pbc$stage4 <- 1*(pbc$stage == 4)
pbc$edema1 <- 1*((pbc$edema == 1)|(pbc$edema == 0.5))
LogNormalReg <- survreg(Surv(time,status) ~ drug + sex + ascites +
                          hepatom + spiders + edema1 + age +
                          bili + chol + albumin + copper + alk + sgot +
                          trig + platelet + prothrombin + stage4,
                        data = pbc,
                        dist = "lognormal")
summary(LogNormalReg)

X <- pbc[,c("drug","sex","ascites","hepatom","spiders","edema1","age","bili","chol","albumin","copper","alk","sgot","trig","platelet","prothrombin","stage4")]

JAGS_AFT = function(){
  # Likelihood
  for(i in 1:276){
    Y[i] ~ dlnorm(mu[i],inv_sigma2)
    mu[i] <- inprod(X[i,],beta) + sigma*W
  }
  # Prior for beta
  for(j in 1:17){
    beta[j] ~ dnorm(0,0.0001)
  }
  # Prior for sigma
  inv_sigma2 ~ dgamma(0.0001,0.0001)
  sigma <- sqrt(1.0/inv_sigma2)
  # Prior for W
  W ~ dnorm(0,1)
}

fit_JAGS_AFT = jags(
  data = list(Y = pbc$time, X = X),
  inits = list(list(inv_sigma2=1,
                    beta=rnorm(17),
                    W=0)),
  parameters.to.save = c("beta"),
  n.chains = 1,
  n.iter = 10000,
  n.burnin = 1000,
  model.file = JAGS_AFT
)
print(fit_JAGS_AFT)