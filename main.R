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
