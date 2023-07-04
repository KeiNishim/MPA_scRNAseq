library(pROC)
#Read Relapse/IFNA/Mono Rate% matrix file

#using GLM
GLM <- glm(outcome ~ IFNA+Mono, family=gaussian, data=df)
ROC <-roc(GLM$y , GLM$fitted.values,ci=T)
auc(ROC)

#Youden Index
plot(ROC,
     identity = T,
     print.thres = T,
     print.thres.best.method="youden",
     legacy.axes = F
)