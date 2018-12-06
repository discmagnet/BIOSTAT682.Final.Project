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
rf1<-randomForest(logTime ~ ., data=pbc2, importance=T)
rf1$importance

library(survival)
library(MASS)
LogNormalReg <- survreg(Surv(pbc2$logTime)~.,data=pbc2, dist = "weibull")
summary(LogNormalReg)

#Cox Model
pbc$SurvObj <- with(pbc, Surv(time, status == 1))
pbc1<-pbc[,-(1:2)]

pbc1$treatment <- as.factor(pbc1$treatment)
pbc1$sex <- as.factor(pbc1$sex)
pbc1$ascites <- as.factor(pbc1$ascites)
pbc1$edema <- as.factor(pbc1$edema)
pbc1$hepatom <- as.factor(pbc1$hepatom)
pbc1$spiders <- as.factor(pbc1$spiders)
pbc1$stage <- as.factor(pbc1$stage)

res.cox1 <- coxph(pbc1$SurvObj~., data =  pbc1)
summary(res.cox1)

# --- Step Function (example with different data)---#
set.seed(625)
# Define a base model - intercept only
HSP.data$Y.level = as.numeric(HSP.data$Y.level)
null.model <- lm(Y.level ~ 1 , data=HSP.data)
# Define the full model - including all predictors
full.model <- lm(Y.level ~ . , data=HSP.data)
HSP_step   <- step(null.model, scope=list(lower=null.model, upper=full.model), direction = 'both', k=2, trace = F)
summary(HSP_step); HSP_step

