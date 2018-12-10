setwd("~/Documents/U Mich Biostats/Biostats 682")

rm(list=ls(all=TRUE))

pbc<-read.csv("https://raw.githubusercontent.com/MLSurvival/ESP/master/ESP_TKDE2016/Dataset/pbc.csv")
pbc$drug <- 1*(pbc$treatment==1)
pbc$female <- 1*(pbc$sex == 1)
pbc$stage4 <- 1*(pbc$stage == 4)
pbc$edema1 <- 1*((pbc$edema == 1)|(pbc$edema == 0.5))

x <- pbc[,c("drug","sex","ascites","hepatom","spiders","edema1","age","bili","chol","albumin","copper","alk","sgot","trig","platelet","prothrombin","stage4")]

temp = scale(x[,7:16])
x_new = cbind(x[,1:6],temp,x[,17])
colnames(x_new)<-colnames(x)

t<-pbc$time
is.na(t)<-pbc$status==0

is.censored<-1-pbc$status
t.cen<-pbc$time+10*pbc$status #+1

tinits1<-pbc$time+500
is.na(tinits1)<-pbc$status==1
#tinits2<-tinits1+5

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
  for(i in 1:n_train) {
    is.censored.train[i] ~ dinterval(t.train[i], t.cen.train[i])
    t.train[i] ~ dlnorm(lambda.train[i], taue)
    lambda.train[i] <- exp(beta0 + inprod(x.train[i,],beta))
  }
  beta0 ~ dnorm(0,1)
  
  for (i in 1:17) {beta[i]  ~ dnorm(0,1)}
  
  for (i in 1:n_test) {
    is.censored.test[i] ~ dinterval(t.pred[i], t.cen.test[i])
    t.pred[i] ~ dlnorm(lambda.pred[i], taue)
    lambda.pred[i] <- exp(beta0 + inprod(x.test[i,],beta))
  }
  taue  ~ dgamma(0.001, 0.001)
}

pbcdata<-list(n_train=as.integer(length(t.train)), 
              n_test =as.integer(length(t.test )), 
              t.train=t.train,
              is.censored.train=is.censored.train, is.censored.test=is.censored.test,
              t.cen.train=t.cen.train, t.cen.test=t.cen.test,
              x.train=as.matrix(x.train), x.test=as.matrix(x.test))

# pbcinits<-list(list(t.train=tinits1.train,
#                       beta0=rnorm(1),
#                       beta =rnorm(17)))

pbcjags_lnorm = jags(
  data = pbcdata,
  #inits = pbcinits,
  parameters.to.save = c("t.pred","beta0", "beta"),
  n.chains = 1,
  n.iter = 50000,
  n.burnin = 1000,
  model.file = surv_model)

print(pbcjags_lnorm)

# DIC 
DIC = pbcjags_lnorm$BUGSoutput$DIC
# pD
pD  = pbcjags_lnorm$BUGSoutput$pD


mcmc_fit <- as.mcmc(pbcjags_lnorm)
#predicted time
t_pred_sample <- mcmc_fit[[1]][,paste("t.pred[",1:nrow(x.test),"]",sep="")]
t_pred=apply(t_pred_sample,2,mean)
t_pred_CI = apply(t_pred_sample,2,quantile,prob=c(0.025,0.975))

t_pred_lcl = t_pred_CI[1,]; t_pred_ucl = t_pred_CI[2,]

#beta
beta_res = mcmc_fit[[1]][,c(paste0("beta[",1:17,"]"), "beta0")]

# remove missing values
t_pred_new <- t_pred[-which(is.na(t.test))]
t.test.new <- as.double(na.omit(t.test))
t_pred_lcl_new <- t_pred_lcl[-which(is.na(t.test))]
t_pred_ucl_new <- t_pred_ucl[-which(is.na(t.test))]
# PMSE
PMSE = mean((t_pred_new-t.test.new)^2)
# coverage
coverage = mean((t.test.new> t_pred_lcl_new)&(t.test.new< t_pred_ucl_new))
