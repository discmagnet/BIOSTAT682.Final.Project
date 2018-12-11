rm(list=ls(all=TRUE))

setwd("~/Downloads/2018_Fall/682/proj")
library(R2jags)

pbc<-read.csv("https://raw.githubusercontent.com/MLSurvival/ESP/master/ESP_TKDE2016/Dataset/pbc.csv")
pbc$drug <- 1*(pbc$treatment==1)
pbc$female <- 1*(pbc$sex == 1)
pbc$stage4 <- 1*(pbc$stage == 4)
pbc$edema1 <- 1*((pbc$edema == 1)|(pbc$edema == 0.5))

x <- pbc[,c("drug","sex","ascites","hepatom","spiders","edema1","age","bili","chol","albumin","copper","alk","sgot","trig","platelet","prothrombin","stage4")]

temp = scale(x[,7:16])
x_new = cbind(x[,1:6],temp,x[,17])
colnames(x_new) <- colnames(x)

t<-(pbc$time)/365
is.na(t)<-pbc$status==0

is.censored<-1-pbc$status
t.cen<-(pbc$time)/365 + pbc$status

tinits1<-(pbc$time)/365 + 5
is.na(tinits1)<-pbc$status==1

set.seed(8102)
train_int <- sample(nrow(pbc),floor(nrow(pbc)*0.8))

is.censored.train = is.censored[train_int]
is.censored.test  = is.censored[-train_int]
t.train = t[train_int]; 
t.test = t[-train_int]
t.cen.train = t.cen[train_int]; 
t.cen.test = t.cen[-train_int]
tinits1.train = tinits1[train_int]; 
tinits1.test = tinits1[-train_int]
x.train = x_new[train_int,]; 
x.test = x_new[-train_int,]


surv_model = function(){
        
        #likelihood
        for(i in 1:n_train) {
                is.censored.train[i] ~ dinterval(t.train[i], t.cen.train[i])
                t.train[i] ~ dweibull(alpha, lambda.train[i])
                log(lambda.train[i]) <- beta0 + inprod(x.train[i,],beta)
        }
        
        # Prior for intercept
        beta0 ~ dnorm(0, 0.0001)
        
        # Prior for beta
        for(j in 1:2) {beta[j] ~ ddexp(0,inv_tau2_demo)}
        beta[3] ~ ddexp(0,inv_tau2_drug)
        for(j in 4:8) {beta[j] ~ ddexp(0,inv_tau2_symp)}
        for(j in 9:16){beta[j] ~ ddexp(0,inv_tau2_sign)}
        beta[17] ~ ddexp(0,inv_tau2_cont)
        
        # Prior for the inverse variance
        inv_tau2_demo ~ dgamma(0.0001,0.0001)
        inv_tau2_drug ~ dgamma(0.0001,0.0001)
        inv_tau2_symp ~ dgamma(0.0001,0.0001)
        inv_tau2_sign ~ dgamma(0.0001,0.0001)
        inv_tau2_cont ~ dgamma(0.0001,0.0001)
        
        #prior for alpha
        alpha ~ dgamma(1.1, 1.1)
        
        #prediction
        for (i in 1:n_test) {
                is.censored.test[i] ~ dinterval(t.pred[i], t.cen.test[i])
                t.pred[i] ~ dweibull(alpha, lambda.pred[i])
                lambda.pred[i] <- exp(beta0 + inprod(x.test[i,],beta))
        }
}


pbcdata<-list(n_train=as.integer(length(t.train)), 
              n_test =as.integer(length(t.test )), 
              t.train=t.train,
              is.censored.train=is.censored.train, is.censored.test=is.censored.test,
              t.cen.train=t.cen.train, t.cen.test=t.cen.test,
              x.train=as.matrix(x.train), x.test=as.matrix(x.test))

pbcjags = jags(
        data = pbcdata,
        #inits = pbcinits,
        parameters.to.save = c("t.pred","beta0", "beta"),
        n.chains = 1,
        n.iter = 50000,
        n.burnin = 1000,
        model.file = surv_model)

mcmc_fit = as.mcmc(pbcjags)
print(pbcjags)
# DIC 
DIC = pbcjags$BUGSoutput$DIC
DIC
# pD
pD  = pbcjags$BUGSoutput$pD
pD

#predicted time
t_pred_sample <- mcmc_fit[[1]][,paste("t.pred[",1:nrow(x.test),"]",sep="")]
t_pred=apply(t_pred_sample,2,mean)
t_pred_CI = apply(t_pred_sample,2,quantile,prob=c(0.025,0.975))

t_pred_lcl = t_pred_CI[1,]; t_pred_ucl = t_pred_CI[2,]

#beta
beta_res = mcmc_fit[[1]][,c(paste0("beta[",1:17,"]"), "beta0")]

# model evaluation
t_pred_new <- t_pred[-which(is.na(t.test))]
t.test.new <- as.double(na.omit(t.test))
t_pred_lcl_new <- t_pred_lcl[-which(is.na(t.test))]
t_pred_ucl_new <- t_pred_ucl[-which(is.na(t.test))]
# PMSE
PMSE = mean((t_pred_new-t.test.new)^2)
PMSE
# coverage
coverage = mean((t.test.new> t_pred_lcl_new)&(t.test.new< t_pred_ucl_new))
coverage

# plots
xyplot(beta_res[,1:3])
densityplot(beta_res)
traceplot(beta_res)
autocorr.plot(beta_res)

