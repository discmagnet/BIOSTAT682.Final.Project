pbc<-read.csv("https://raw.githubusercontent.com/MLSurvival/ESP/master/ESP_TKDE2016/Dataset/pbc.csv")

library("glmnet", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(randomForest)
pbc$logTime<-log(pbc$time)
pbc2<-subset(pbc, pbc$status==1)
pbc2 <- pbc2[,-(1:2)]
X <- as.matrix(pbc2[,-c(1,2,20)])
Y <- as.matrix(pbc2[,20])


pbc2$treatment <- as.factor(pbc2$treatment)
pbc2$sex <- as.factor(pbc2$sex)
pbc2$ascites <- as.factor(pbc2$ascites)
pbc2$edema <- as.factor(pbc2$edema)
pbc2$hepatom <- as.factor(pbc2$hepatom)
pbc2$spiders <- as.factor(pbc2$spiders)
pbc2$stage <- as.factor(pbc2$stage)

#BLASSO<-blasso(X, Y, T = 1000, thin = 1, RJ = TRUE, M = NULL,
#       beta = NULL, lambda2 = 1, s2 = var(Y-mean(Y)),
#      icept = TRUE, 
#       normalize = TRUE, verb = 1)


#randomForest::randomForest(x = pbc2[,-c(1,2,20)],y=as.factor(pbc2[,20]))
#rf_res <- randomForest::randomForest(x = pbc2[,-c(1,2,20)],y=as.factor(pbc2[,20]))
#predict(rf_res)
<<<<<<< HEAD
rf1<-randomForest::randomForest(logTime ~ ., data=pbc2, importance=T)
=======
rf1<-randomForest(logTime ~ ., data=pbc2, importance=T)
rf1$importance

library(survival)
library(MASS)
LogNormalReg <- survreg(Surv(pbc2$logTime)~.,data=pbc2, dist = "weibull")
summary(LogNormalReg)

# --- Cox model ---
pbc$time <- with(pbc, Surv(time, status == 1))

pbc$treatment <- as.factor(pbc$treatment)
pbc$sex <- as.factor(pbc$sex)
pbc$edema <- as.factor(pbc$edema)
pbc$ascites <- as.factor(pbc$ascites)
pbc$hepatom <- as.factor(pbc$hepatom)
pbc$spiders <- as.factor(pbc$spiders)
pbc$stage <- as.factor(pbc$stage)

pbc.cox <- pbc[,-2]
full.cox <- coxph(time ~ ., data =  pbc.cox)
summary(full.cox)

null.cox <- coxph(time ~ 1, data =  pbc.cox)
summary(null.cox)
# res.zph1 <- cox.zph(res.cox1)
# plot(res.zph1)

cox_step <- step(null.cox, scope=list(lower=null.cox, upper=full.cox), direction = 'both', k=2, trace = F)
summary(cox_step)

>>>>>>> 2524aee0516618cc60eeec855409beb014571962
