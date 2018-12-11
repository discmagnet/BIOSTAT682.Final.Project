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

rf1<-randomForest::randomForest(logTime ~ ., data=pbc2, importance=T)

rf1<-randomForest(logTime ~ ., data=pbc2, importance=T)
rf1$importance

library(survival)
library(MASS)
LogNormalReg <- survreg(Surv(pbc2$logTime)~.,data=pbc2, dist = "weibull")
summary(LogNormalReg)

# --- AFT model ---


pbc$treatment <- as.factor(pbc$treatment)
pbc$edema <- as.factor(pbc$edema)
pbc$ascites <- as.factor(pbc$ascites)
pbc$hepatom <- as.factor(pbc$hepatom)
pbc$spiders <- as.factor(pbc$spiders)
pbc$stage <- as.factor(pbc$stage)
pbc$female <- 1*(pbc$sex == 1)
pbc$stage4 <- 1*(pbc$stage == 4)
pbc$edema1 <- 1*((pbc$edema == 1)|(pbc$edema == 0.5))

full.aft <- survreg(Surv(time,status) ~ treatment + age + female + ascites + hepatom +
                      spiders + edema1 + bili + chol + albumin + copper + alk + sgot +
                      trig + platelet + prothrombin + stage4,
                    data = pbc,
                    dist = "lognormal")
#summary(full.aft)

null.aft <- survreg(Surv(time,status) ~ 1,
                    data = pbc,
                    dist = "lognormal")
#summary(null.aft)
# res.zph1 <- cox.zph(res.cox1)
# plot(res.zph1)

aft_step <- step(null.aft, scope=list(lower=null.aft, upper=full.aft), direction = 'both', k=2, trace = F)
summary(aft_step)

# Look at the baseline distribution of survival times
library(ggplot2)
plot01 <- ggplot(data = pbc, aes(x=time)) + 
  geom_histogram(show.legend = FALSE, aes(y = ..density..), 
    bins = 20,colour = "black",fill = "white") +
  geom_density(colour = "#4271AE") +
  labs(x = "Time (in days)", y = "Freq") +
  ggtitle("Distribution of Survival Times") +
  theme(plot.title = element_text(hjust = 0.5)) 
plot01

plot02 <- ggplot(data = pbc, aes(x=logTime)) + 
  geom_histogram(show.legend = FALSE, aes(y = ..density..), 
                 bins = 20,colour = "black",fill = "white") +
  geom_density(colour = "#4271AE") +
  labs(x = "Time (in log scale)", y = "Freq") +
  ggtitle("Distribution of Log(Ti)") +
  theme(plot.title = element_text(hjust = 0.5)) 
plot02

pbc$Years <- pbc$time/365
pbc$logYears <- log(pbc$time/365)

plot03 <- ggplot(data = pbc, aes(x=Years)) + 
  geom_histogram(show.legend = FALSE, aes(y = ..density..), 
                 bins = 20,colour = "black",fill = "white") +
  geom_density(colour = "#4271AE") +
  labs(x = "Time (in years)", y = "Freq") +
  ggtitle("Distribution of Survival Times") +
  theme(plot.title = element_text(hjust = 0.5)) 
plot03

plot04 <- ggplot(data = pbc, aes(x=logYears)) + 
  geom_histogram(show.legend = FALSE, aes(y = ..density..), 
                 bins = 20,colour = "black",fill = "white") +
  geom_density(colour = "#4271AE") +
  labs(x = "Time (in log scale)", y = "Freq") +
  ggtitle("Distribution of Log(Ti)") +
  theme(plot.title = element_text(hjust = 0.5)) 
plot04

# ploting simple KM curves

pbc$SurvObj <- with(pbc,Surv(time,status == 1))
km.as.one <- survfit(SurvObj ~ 1, data = pbc, conf.type = "log-log")
plot(km.as.one)

km.by.sex <- survfit(SurvObj ~ 1, data = pbc, conf.type = "log-log")