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
                t.train[i] ~ dexp(lambda.train[i])
                log(lambda.train[i]) <- beta0 + inprod(x.train[i,],beta)
        }
        beta0 ~ dnorm(0,1)
        
        for (i in 1:17) {beta[i]  ~ dnorm(0,1)}
        
        for (i in 1:n_test) {
                is.censored.test[i] ~ dinterval(t.pred[i], t.cen.test[i])
                t.pred[i] ~ dexp(lambda.pred[i])
                lambda.pred[i] <- exp(beta0 + inprod(x.test[i,],beta))
        }
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

pbcjags = jags(
        data = pbcdata,
        #inits = pbcinits,
        parameters.to.save = c("t.pred","beta0", "beta"),
        n.chains = 1,
        n.iter = 50000,
        n.burnin = 1000,
        model.file = surv_model)

print(pbcjags)
