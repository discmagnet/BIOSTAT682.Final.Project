rm(list=ls(all=TRUE))

setwd("~/Downloads/2018_Fall/682/proj")
library(R2jags)

pbc<-read.csv("https://raw.githubusercontent.com/MLSurvival/ESP/master/ESP_TKDE2016/Dataset/pbc.csv")
pbc$drug <- 1*(pbc$treatment==1)
pbc$female <- 1*(pbc$sex == 1)
pbc$stage4 <- 1*(pbc$stage == 4)
pbc$edema1 <- 1*((pbc$edema == 1)|(pbc$edema == 0.5))

Y = pbc$time
X = pbc
X$time = rep(1, nrow(X))
colnames(X)[1] = "intercept"
X <- X[,c("age","sex","drug","ascites","hepatom","spiders","edema1","stage4","bili","chol","albumin","copper","alk","sgot","trig","platelet","prothrombin")]
# Demographics:   age sex 
# Treatment: drug
# symptom:   status ascites hepatom spiders edema1 stage4
# Vital signs： bili chol albumin copper alk sgot trig platelet 
# continuous:   protime

#split data into train and test
set.seed(8102)
train_int <- sample(nrow(pbc),floor(nrow(pbc)*0.8))

Y_train = Y[train_int ]
Y_test  = Y[-train_int]
X_train = X[train_int, ]
X_test  = X[-train_int,]

# explore data
par(mfcol=c(1,2))
plot(X_train$bili,Y_train,col="blue",xlab="serum bilirubin in mg/dl",
     ylab="log(time)",cex.lab=1.5,cex.axis=1.5,
     ylim=c(0,5000),main="Training Data")

plot(X_test$bili,Y_test,col="red",xlab="serum bilirubin in mg/dl",
     ylab="log(time)",cex.lab=1.5,cex.axis=1.5,
     ylim=c(0,5000),main="Test Data")


# --- JAGS_AFT function: BLASSO --- #
JAGS_BLASSO = function(Y_train,X_train,X_test,n.iter=10000,n.burnin=1000){
    
    JAGS_AFT = function(){
            # Likelihood
            for(i in 1:n_train){
                    Y_train[i] ~ dlnorm(mu[i],inv_sigma2)
                    mu[i] <- beta0 + inprod(X_train[i,],beta) + sigma*W
            }
            # Prior for intercept
            beta0 ~ dnorm(0, 0.0001)
            
            # Prior for beta
            for(j in 1:2) {beta[j] ~ ddexp(0,inv_tau2_demo)}
                           beta[3] ~ ddexp(0,inv_tau2_drug)
            for(j in 4:8) {beta[j] ~ ddexp(0,inv_tau2_symp)}
            for(j in 9:16){beta[j] ~ ddexp(0,inv_tau2_sign)}
                           beta[p] ~ ddexp(0,inv_tau2_cont)
        
            # Prior for the inverse variance
            inv_tau2_demo ~ dgamma(0.0001,0.0001)
            inv_tau2_drug ~ dgamma(0.0001,0.0001)
            inv_tau2_symp ~ dgamma(0.0001,0.0001)
            inv_tau2_sign ~ dgamma(0.0001,0.0001)
            inv_tau2_cont ~ dgamma(0.0001,0.0001)
        
            inv_sigma2 ~ dgamma(0.0001, 0.0001)
            sigma <- sqrt(1.0/inv_sigma2)
                       
            #prior for W
            W ~ dnorm(0,1)
        
            #prediction
            for (i in 1:n_pred) {
                    Y_pred[i] ~ dlnorm(mu_pred[i],inv_sigma2)
                    mu_pred[i] <- beta0 + inprod(X_test[i,],beta) + sigma*W
            }
    }

    AFT.data = list(Y_train=Y_train,
                    X_train=X_train,
                    X_test =X_test ,
                    n_train=as.integer(nrow(X_train)),
                    n_pred =as.integer(nrow(X_test)),
                    p=ncol(X_train)
    )

    fit_JAGS_BLASSO = jags(
                       data=AFT.data,
                       inits=list(list(beta = rnorm(17),
                                       beta0 = rnorm(1),
                                       inv_sigma2 = 1,
                                       inv_tau2_demo = 1,
                                       inv_tau2_drug = 1,
                                       inv_tau2_symp = 1,
                                       inv_tau2_sign = 1,
                                       inv_tau2_cont = 1,
                                       W=0)),
                       parameters.to.save = c("Y_pred", "beta0", "beta"),
                       n.chains=1,
                       n.iter=10000,
                       n.burnin=1000,
                       model.file=JAGS_AFT)
    
        
    mcmc_fit = as.mcmc(fit_JAGS_BLASSO)
    #predicted Y
    Y_pred_sample = mcmc_fit[[1]][,paste("Y_pred[",1:nrow(X_test),"]",sep="")]
    Y_pred=apply(Y_pred_sample,2,mean)
    Y_pred_CI = apply(Y_pred_sample,2,quantile,prob=c(0.025,0.975))
    #beta
    beta_res = mcmc_fit[[1]][,c(paste0("beta[",1:17,"]"), "beta0")]
    return(list(
                Y_pred=Y_pred,
                beta_res=beta_res,
                Y_pred_lcl = Y_pred_CI[1,],
                Y_pred_ucl = Y_pred_CI[2,],
                DIC = fit_JAGS_BLASSO$BUGSoutput$DIC,
                pD  = fit_JAGS_BLASSO$BUGSoutput$pD,
                fit.JAGS=fit_JAGS_BLASSO))
}


plot_res = function(res,main=""){
        plot(X_test$age,res$Y_pred_ucl,type="l",
             ylim=c(0,5000),xlab="age/day",ylab="log(Time)",
             cex.lab=1.5,cex.axis=1.5,main=main)
        lines(X_test$age,res$Y_pred_lcl)
        lines(X_test$age,res$Y_pred,col="blue")
        points(X_test$age, Y_test,col="red")
}

summary_res = function(res,Y_test){
        PMSE = mean((res$Y_pred-Y_test)^2)
        coverage = mean((Y_test>res$Y_pred_lcl)&(Y_test<res$Y_pred_ucl))
        return(c(PMSE=PMSE,coverage=coverage))
}

res_aft = JAGS_BLASSO( Y_train=Y_train,
                       X_train=X_train,
                       X_test=X_test)

# check result
plot_res(res_aft,main=sprintf("BLASSO, DIC = %.2f",res_aft$DIC))
summary_res(res_aft, Y_test)
# plots
xyplot(res_aft$beta_res[,1:3])
densityplot(mcmc)
traceplot(res_aft$beta_res)
autocorr.plot(res_aft$beta_res)
