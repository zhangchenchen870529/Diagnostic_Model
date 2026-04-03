#####加载packages
source("./rfcv1.R")
library(openxlsx)
library(readr)
library(data.table)
library(tibble)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(bitops)
library(Rcpp)
library(foreign)
library(VIM)
library(seqinr)
library(edgeR)
library(limma)
library(DESeq2)
library(snowfall)
library(doParallel)
library(future.apply)
library(caret)
library(randomForest)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(mboost)
library(e1071)
library(MASS)
library(xgboost)
library(gbm)
library(kernlab)
library(neuralnet)
library(NeuralNetTools)
library(BART)
library(obliqueRSF)
library(aorsf)
library(party)
library(partykit)
library(Metrics)
library(DALEX)
library(pROC)
library(rpart)
library(rpart.plot)
library(rattle)
library(survival)
library(survminer)
library(plsRcox)
library(superpc)
library(CoxBoost)
library(survivalsvm)
library(rms)
library(pec)
library(ggDCA)
library(timeROC)
library(regplot)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(ggbreak)
library(ggsignif)
library(ggpol)
library(sparkline)
library(visNetwork)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(venn)
library(VennDiagram)
library(UpSetR)
library(beepr)

argv=commandArgs(T)

####训练集及测试集
data=read.table(argv[1],header = T,sep = "\t",check.names = F,row.names = 1)
data=as.matrix(data)

train=read.table(argv[2],header = T,sep = "\t",check.names = F)
id00=train[train[,2]==0,]
id01=train[train[,2]==1,]
id0=rbind(id00,id01)
grouptrain=id0[,2]
trainexp=t(data[,pmatch(id0[,1],colnames(data)),drop=F])

test1=read.table(argv[3],header = T,sep = "\t",check.names = F)
id10=test1[test1[,2]==0,]
id11=test1[test1[,2]==1,]
id1=rbind(id10,id11)
grouptest1=id1[,2]
test1exp=t(data[,pmatch(id1[,1],colnames(data)),drop=F])

test2=read.table(argv[4],header = T,sep = "\t",check.names = F)
id20=test2[test2[,2]==0,]
id21=test2[test2[,2]==1,]
id2=rbind(id20,id21)
grouptest2=id2[,2]
test2exp=t(data[,pmatch(id2[,1],colnames(data)),drop=F])

test3=read.table(argv[5],header = T,sep = "\t",check.names = F)
id30=test3[test3[,2]==0,]
id31=test3[test3[,2]==1,]
id3=rbind(id30,id31)
grouptest3=id3[,2]
test3exp=t(data[,pmatch(id3[,1],colnames(data)),drop=F])

grouptrain
trainexp[1:5,1:5]

train=data.frame(group=grouptrain,trainexp)
test1=data.frame(group=grouptest1,test1exp)
test2=data.frame(group=grouptest2,test2exp)
test3=data.frame(group=grouptest3,test3exp)

datalist=list(Train=train,Test1=test1,Test2=test2,Test3=test3)

result=data.frame()
modelgenes=list()
dir.create("ROC")

#####################################################
######1.Enet ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                   y = grouptrain,
                   family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
             y = grouptrain,
             family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  modelname=paste0('Enet','[alpha=',alpha,']')
  modelgenes[[modelname]]=Enetgenes
  
  rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
  AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]})) %>% rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
  pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
        text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1) 
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######1.2Enet+glm ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)
  fit <- step(glm(formula = ifelse(grouptrain=="0",0,1) ~ .,
                  family = "binomial", 
                  data = as.data.frame(trainexp2)),trace = 0)
  fit$subFeature = colnames(trainexp2)
  glmgenes=names(coef(fit))[2:length(names(coef(fit)))]
  modelname=paste0('Enet','[alpha=',alpha,']+glm')
  modelgenes[[modelname]]=glmgenes
  rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = 'response', as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
  AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]})) %>% rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
  pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+glm',"_AUC_result"),sep="\t",quote=F,col.names=NA)
#####################################################
######1.3Enet+svm ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = as.factor(grouptrain),
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = as.factor(grouptrain),
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

cl <- makeCluster(16)
registerDoParallel(cl)
fit=rfe(x=as.data.frame(trainexp2),
        y=as.factor(grouptrain),
        rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
        methods="svmRadial")
stopCluster(cl)
modelname=paste0('Enet','[alpha=',alpha,']+svm')
modelgenes[[modelname]]=fit[["optVariables"]]

rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
        text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1) 
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+svm',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######1.4Enet+XGBoost ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]

  trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

dtrain <- xgb.DMatrix(data = as.matrix(trainexp2), label=ifelse(grouptrain=="0",0,1))
xgb_params <- list(objective = "binary:logistic", eval_metric = "auc",base_score = 0.5,max_depth = 3,eta = 0.1,subsample = 0.8,colsample_bytree = 0.8)
initial_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = 100,verbose = 0)

feature_importance <- xgb.importance(feature_names = colnames(trainexp2), model = initial_model)
selected_features <- feature_importance$Feature  # 按重要性降序排列

selected_features
modelname=paste0('Enet','[alpha=',alpha,']+XGBoost')
modelgenes[[modelname]] <- selected_features

trainexp_selected <- trainexp2[, selected_features, drop = FALSE]
dtrain_selected <- xgb.DMatrix(data = trainexp_selected, label = grouptrain)

final_model <- xgb.train(params = xgb_params,data = dtrain_selected,nrounds = 100,verbose = 0)

rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(final_model, type = "response", as.data.frame(x[,selected_features,drop=F])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS!="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS!="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+XGBoost',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######1.5Enet+RF  ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = grouptrain,
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = grouptrain,
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
   trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)
  

set.seed(as.numeric(argv[6]))
  train.cv <- replicate(5, rfcv1(trainexp2, as.factor(grouptrain), cv.fold = 5, step = 0.9), simplify = F)
  error.cv <- sapply(train.cv, "[[", "error.cv")
  error.cv.rm <- rowMeans(error.cv)
  id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  marker.t <- sort(marker.t, d = T)
  names(marker.t) <- colnames(trainexp2)[as.numeric(names(marker.t))]
  rfGenes <- names(marker.t)[1:marker.num]
  modelname=paste0('Enet','[alpha=',alpha,']+RF')
  modelgenes[[modelname]]=rfGenes

  rf2 <- randomForest(as.factor(grouptrain)~.,data=data.frame(trainexp2[,rfGenes,drop=F]), importance = T)
  rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(rf2, type = "prob", as.data.frame(x[,2:ncol(x)])))[,2])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
  AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
  pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+RF',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######1.6Enet+rpart ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = as.factor(grouptrain),
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = as.factor(grouptrain),
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)
mod1<-rpart(as.factor(grouptrain)~.,data = as.data.frame(trainexp2),method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>% arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('Enet','[alpha=',alpha,']+rpart')
modelgenes[[modelname]]=JCSGenes
print (JCSGenes)
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()
 AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+rpart',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######1.7Enet+GBM ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = as.factor(grouptrain),
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = as.factor(grouptrain),
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0.1]
modelname=paste0('Enet','[alpha=',alpha,']+GBM')
modelgenes[[modelname]]=GBMgenes

rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
 pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+GBM',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######1.8Enet+Ridge ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = as.factor(grouptrain),
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = as.factor(grouptrain),
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('Enet','[alpha=',alpha,']+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+ridge',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######1.9Enet+Lasso ######
#####################################################
set.seed(as.numeric(argv[6]))
for (alpha in seq(0.1,0.9,0.1)) {
  cv.fit = cv.glmnet(x = trainexp,
                     y = as.factor(grouptrain),
                     family = "binomial", alpha = alpha, nfolds = 100)
  fit = glmnet(x = trainexp,
               y = as.factor(grouptrain),
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(trainexp)
  Enetgenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
  Enetgenes=Enetgenes[2:length(Enetgenes)]
  trainexp2=as.data.frame(trainexp[,Enetgenes])
  test1exp2=as.data.frame(test1exp[,Enetgenes])
  test2exp2=as.data.frame(test2exp[,Enetgenes])
  test3exp2=as.data.frame(test3exp[,Enetgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('Enet','[alpha=',alpha,']+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
}
write.table(result,paste0("./ROC/",'Enet+lasso',"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######2.glm   ######
#####################################################
fit <- step(glm(formula = ifelse(grouptrain=="0",0,1) ~ .,
                family = "binomial", 
                data = as.data.frame(trainexp)),trace = 0)
fit$subFeature = colnames(trainexp)
glmgenes=names(coef(fit))[2:length(names(coef(fit)))]
modelname=paste0('glm')
modelgenes[[modelname]]=glmgenes
rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = 'response', as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()
AUCs$Model <- modelname
result <- rbind(result,AUCs)
write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

  trainexp2=as.data.frame(trainexp[,glmgenes])
  test1exp2=as.data.frame(test1exp[,glmgenes])
  test2exp2=as.data.frame(test2exp[,glmgenes])
  test3exp2=as.data.frame(test3exp[,glmgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

#####################################################
######2.2.glm+SVM   ######
#####################################################
cl <- makeCluster(16)
registerDoParallel(cl)
fit=rfe(x=trainexp2,
        y=as.factor(grouptrain),
        rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
        methods="svmRadial")
stopCluster(cl)
modelname=paste0('glm','+svm')
modelgenes[[modelname]]=fit[["optVariables"]]
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######2.3.glm+XGBoost   ######
#####################################################
# 训练XGBoost模型
dtrain <- xgb.DMatrix(data = as.matrix(trainexp2), label=ifelse(grouptrain=="0",0,1))
xgb_params <- list(objective = "binary:logistic", eval_metric = "auc",base_score = 0.5,max_depth = 3,eta = 0.1,subsample = 0.8,colsample_bytree = 0.8)
initial_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = 100,verbose = 0)

# 获取特征重要性并排序
feature_importance <- xgb.importance(feature_names = colnames(trainexp2), model = initial_model)
selected_features <- feature_importance$Feature  # 按重要性降序排列

selected_features
modelname=paste0('glm','+XGBoost')
modelgenes[[modelname]] <- selected_features

# 4. 使用最佳特征子集训练最终模型
trainexp_selected <- trainexp2[, selected_features, drop = FALSE]
#x_test_selected <- x_test[, selected_features, drop = FALSE]
dtrain_selected <- xgb.DMatrix(data = trainexp_selected, label = grouptrain)

# 训练最终模型
final_model <- xgb.train(params = xgb_params,data = dtrain_selected,nrounds = 100,verbose = 0)

rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(final_model, type = "response", as.data.frame(x[,selected_features,drop=F])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
#####################################################
######2.4.glm+RF  ######
#####################################################
set.seed(0)
  train.cv <- replicate(5, rfcv1(trainexp2, as.factor(grouptrain), cv.fold = 5, step = 0.9), simplify = F)
  error.cv <- sapply(train.cv, "[[", "error.cv")
  error.cv.rm <- rowMeans(error.cv)
  id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
  
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
 
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))

  marker.t <- sort(marker.t, d = T)
 
  names(marker.t) <- colnames(trainexp2)[as.numeric(names(marker.t))]
  rfGenes <- names(marker.t)[1:marker.num]
  modelname=paste0('glm','+RF')
  modelgenes[[modelname]]=rfGenes

  rf2 <- randomForest(as.factor(grouptrain)~.,data=data.frame(trainexp2[,rfGenes,drop=F]), importance = T) 
  rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(rf2, type = "prob", as.data.frame(x[,2:ncol(x)])))[,2])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
 pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######2.5.glm+rpart  ######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('glm','+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(mod1,as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######2.6.glm+GBM ######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0.1]
modelname=paste0('glm','+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######2.7.glm+Ridge ######
#####################################################
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('glm','+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######2.8.glm+Lasso ######
#####################################################
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
lassogenes=lassogenes[2:length(lassogenes)]
modelname=paste0('glm','+Lasso')
modelgenes[[modelname]]=lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######3.SVM ######
#####################################################
data <- as.data.frame(trainexp)
cl <- makeCluster(16)
registerDoParallel(cl)
fit=rfe(x=data,
        y=as.factor(grouptrain),
        rfeControl=rfeControl(functions=caretFuncs, method="cv"),
        methods="svmRadial")
fit
stopCluster(cl)
modelname=paste0("SVM")
svmgenes=fit[["optVariables"]]
modelgenes[[modelname]]=svmgenes
svmgenes
rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)
write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

  trainexp2=as.data.frame(trainexp[,svmgenes])
  test1exp2=as.data.frame(test1exp[,svmgenes])
  test2exp2=as.data.frame(test2exp[,svmgenes])
  test3exp2=as.data.frame(test3exp[,svmgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

#####################################################
######3.2.SVM+XGBoost ######
#####################################################
dtrain <- xgb.DMatrix(data = as.matrix(trainexp2), label=ifelse(grouptrain=="0",0,1))
xgb_params <- list(objective = "binary:logistic", eval_metric = "auc",base_score = 0.5,max_depth = 3,eta = 0.1,subsample = 0.8,colsample_bytree = 0.8)
initial_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = 100,verbose = 0)

feature_importance <- xgb.importance(feature_names = colnames(trainexp2), model = initial_model)
selected_features <- feature_importance$Feature  # 按重要性降序排列

selected_features
modelname=paste0('SVM+XGBoost')
modelgenes[[modelname]] <- selected_features

trainexp_selected <- trainexp2[, selected_features, drop = FALSE]
dtrain_selected <- xgb.DMatrix(data = trainexp_selected, label = grouptrain)

final_model <- xgb.train(params = xgb_params,data = dtrain_selected,nrounds = 100,verbose = 0)

rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(final_model, type = "response", as.data.frame(x[,selected_features,drop=F])))[,1])})

  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######3.3.SVM+RF ######
#####################################################
set.seed(as.numeric(argv[6]))
  train.cv <- replicate(5, rfcv1(trainexp2, as.factor(grouptrain), cv.fold = 5, step = 0.9), simplify = F)
  error.cv <- sapply(train.cv, "[[", "error.cv")
  error.cv.rm <- rowMeans(error.cv)
  id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
  
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
 
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))

  marker.t <- sort(marker.t, d = T)
 
  names(marker.t) <- colnames(trainexp2)[as.numeric(names(marker.t))]
  rfGenes <- names(marker.t)[1:marker.num]
  modelname=paste0('SVM','+RF')
  modelgenes[[modelname]]=rfGenes

  rf2 <- randomForest(as.factor(grouptrain)~.,data=data.frame(trainexp2[,rfGenes,drop=F]), importance = T) 
  rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(rf2, type = "prob", as.data.frame(x[,2:ncol(x)])))[,2])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######3.4.SVM+rpart ######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('SVM+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######3.5.SVM+GBM ######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0.1]
modelname=paste0('SVM+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######3.6.SVM+Ridge ######
#####################################################
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('SVM+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######3.6.SVM+Lasso ######
#####################################################
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
lassogenes=lassogenes[2:length(lassogenes)]
modelname=paste0('SVM+Lasso')
modelgenes[[modelname]]=lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######4.XGBoost ######
#####################################################
set.seed(as.numeric(argv[6]))
dtrain <- xgb.DMatrix(data = as.matrix(trainexp), label=ifelse(grouptrain=="0",0,1))
xgb_params <- list(objective = "binary:logistic", eval_metric = "auc",base_score = 0.5,max_depth = 3,eta = 0.1,subsample = 0.8,colsample_bytree = 0.8)
initial_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = 100,verbose = 0)

feature_importance <- xgb.importance(feature_names = colnames(trainexp), model = initial_model)
selected_features <- feature_importance$Feature  # 按重要性降序排列

selected_features
modelname=paste0('XGBoost')
modelgenes[[modelname]] <- selected_features

trainexp_selected <- trainexp[, selected_features, drop = FALSE]
dtrain_selected <- xgb.DMatrix(data = trainexp_selected, label = grouptrain)

final_model <- xgb.train(params = xgb_params,data = dtrain_selected,nrounds = 100,verbose = 0)

rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(final_model, type = "response", as.data.frame(x[,selected_features,drop=F])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

  trainexp2=as.data.frame(trainexp[,selected_features])
  test1exp2=as.data.frame(test1exp[,selected_features])
  test2exp2=as.data.frame(test2exp[,selected_features])
  test3exp2=as.data.frame(test3exp[,selected_features])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

#####################################################
######4.2.XGBoost+RF ######
#####################################################
set.seed(as.numeric(argv[6]))
  train.cv <- replicate(5, rfcv1(trainexp2, as.factor(grouptrain), cv.fold = 5, step = 0.9), simplify = F)
  error.cv <- sapply(train.cv, "[[", "error.cv")
  error.cv.rm <- rowMeans(error.cv)
  id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
  
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
 
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))

  marker.t <- sort(marker.t, d = T)
 
  names(marker.t) <- colnames(trainexp2)[as.numeric(names(marker.t))]
  rfGenes <- names(marker.t)[1:marker.num]
  modelname=paste0('XGBoost+RF')
  modelgenes[[modelname]]=rfGenes

  rf2 <- randomForest(as.factor(grouptrain)~.,data=data.frame(trainexp2[,rfGenes,drop=F]), importance = T) 
  rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(rf2, type = "prob", as.data.frame(x[,2:ncol(x)])))[,2])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######4.3.XGBoost+rpart ######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('XGBoost+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######4.4.XGBoost+GBM ######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0.1]
modelname=paste0('XGBoost+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######4.5.XGBoost+Ridge ######
#####################################################
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('XGBoost+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######4.6.XGBoost+Lasso ######
#####################################################
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)


fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('XGBoost+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######5.RF ######
#####################################################
set.seed(as.numeric(argv[6]))
  train.cv <- replicate(5, rfcv1(trainexp, as.factor(grouptrain), cv.fold = 5, step = 0.9), simplify = F)
  error.cv <- sapply(train.cv, "[[", "error.cv")
  error.cv.rm <- rowMeans(error.cv)
  id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))

  marker.t <- sort(marker.t, d = T)
 
  names(marker.t) <- colnames(trainexp)[as.numeric(names(marker.t))]
  rfGenes <- names(marker.t)[1:marker.num]
  modelname=paste0('RF')
  modelgenes[[modelname]]=rfGenes

  rf2 <- randomForest(as.factor(grouptrain)~.,data=data.frame(trainexp[,rfGenes,drop=F]), importance = T) 
  rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(rf2, type = "prob", as.data.frame(x[,2:ncol(x)])))[,2])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

imp=rf2$importance
write.table(result,paste0("./ROC/",modelname,"_feature_importance.txt"),sep="\t",quote=F,col.names=NA)


  trainexp2=as.data.frame(trainexp[,rfGenes])
  test1exp2=as.data.frame(test1exp[,rfGenes])
  test2exp2=as.data.frame(test2exp[,rfGenes])
  test3exp2=as.data.frame(test3exp[,rfGenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

#####################################################
######5.2.RF+rpart ######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = trainexp2,method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('RF+rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######5.3.RF+GBM ######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0.1]
modelname=paste0('RF+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######5.4.RF+Ridge ######
#####################################################
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('RF+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######5.5.RF+Lasso ######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('RF+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######6.rpart ######
#####################################################
mod1<-rpart(as.factor(grouptrain)~.,data = as.data.frame(trainexp),method = "class")
importances <- varImp(mod1)
jcsimportances=as.matrix(importances %>%
                           arrange(desc(Overall)))
JCSGenes=jcsimportances[order(jcsimportances[,"Overall"], decreasing = TRUE),]
JCSGenes=names(JCSGenes[JCSGenes>0])
modelname=paste0('rpart')
modelgenes[[modelname]]=JCSGenes
rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(mod1, as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

  trainexp2=as.data.frame(trainexp[,JCSGenes])
  test1exp2=as.data.frame(test1exp[,JCSGenes])
  test2exp2=as.data.frame(test2exp[,JCSGenes])
  test3exp2=as.data.frame(test3exp[,JCSGenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

#####################################################
######6.2.rpart+GBM######
#####################################################
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp2),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp2)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0.1]
modelname=paste0('rpart+GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######6.3.rpart+Ridge######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('rpart+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######6.4.rpart+Lasso######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('rpart+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######7.GBM######
#####################################################
set.seed(as.numeric(argv[6]))
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = ifelse(grouptrain=="0",0,1) ~ .,
           data = as.data.frame(trainexp),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(trainexp)
GBMgenes=rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0.1]
modelname=paste0('GBM')
modelgenes[[modelname]]=GBMgenes
rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=as.data.frame(predict(fit, type = "response", as.data.frame(x[,2:ncol(x)])))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

  trainexp2=as.data.frame(trainexp[,GBMgenes])
  test1exp2=as.data.frame(test1exp[,GBMgenes])
  test2exp2=as.data.frame(test2exp[,GBMgenes])
  test3exp2=as.data.frame(test3exp[,GBMgenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

#####################################################
######7.2.GBM+Ridge######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('GBM+Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######7.3.GBM+Lasso######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('GBM+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######8.Ridge######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = trainexp,
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 0, nfolds = 1000)
fit = glmnet(x = trainexp,
             y = as.factor(grouptrain),
             family = "binomial", alpha = 0, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp)
Ridgegenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Ridgegenes=Ridgegenes[2:length(Ridgegenes)]
modelname=paste0('Ridge')
modelgenes[[modelname]]=Ridgegenes
rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

  trainexp2=as.data.frame(trainexp[,Ridgegenes])
  test1exp2=as.data.frame(test1exp[,Ridgegenes])
  test2exp2=as.data.frame(test2exp[,Ridgegenes])
  test3exp2=as.data.frame(test3exp[,Ridgegenes])
  train2=data.frame(group=grouptrain,trainexp2)
  test12=data.frame(group=grouptest1,test1exp2)
  test22=data.frame(group=grouptest2,test2exp2)
  test32=data.frame(group=grouptest3,test3exp2)
  datalist2=list(Train=train2,Test1=test12,Test2=test22,Test3=test32)

#####################################################
######8.2.Ridge+Lasso######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = as.matrix(trainexp2),
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = as.matrix(trainexp2),
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp2)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('Ridge+Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist2,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)

#####################################################
######9.Lasso######
#####################################################
set.seed(as.numeric(argv[6]))
cv.fit = cv.glmnet(x = trainexp,
                   y = as.factor(grouptrain),
                   family = "binomial", alpha = 1, nfolds = 1000)
fit = glmnet(x = trainexp,
             y = as.factor(grouptrain),
             family = "binomial", alpha = 1, lambda = cv.fit$lambda.min)
fit$subFeature = colnames(trainexp)
Lassogenes=rownames(coef(fit))[which(coef(fit)[, 1]!=0)]
Lassogenes=Lassogenes[2:length(Lassogenes)]
modelname=paste0('Lasso')
modelgenes[[modelname]]=Lassogenes
rs <- lapply(datalist,function(x){cbind(x[,1,drop=F],RS=predict(fit, type = 'response', as.matrix(x[,2:ncol(x)]))[,1])})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test1"]]=rs[["Test1"]][rs[["Test1"]]$RS!="-Inf",]
  rs[["Test2"]]=rs[["Test2"]][rs[["Test2"]]$RS!="-Inf",]
  rs[["Test3"]]=rs[["Test3"]][rs[["Test3"]]$RS!="-Inf",]
  write.table(rs[["Train"]],paste0("./ROC/",modelname,"_Train_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test1"]],paste0("./ROC/",modelname,"_Test1_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test2"]],paste0("./ROC/",modelname,"_Test2_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(rs[["Test3"]],paste0("./ROC/",modelname,"_Test3_predict_result"),quote=F,sep="\t",col.names=NA)
  write.table(modelgenes[[modelname]],paste0("./ROC/",modelname,"_predict_feature"),quote=F,sep="\t",col.names=NA)

set.seed(as.numeric(argv[6]))
AUCs <- data.frame(auc=sapply(rs,function(x){as.numeric(ci.auc(roc(x$group, as.numeric(x$RS)), method="bootstrap"))[2]}))%>%
  rownames_to_column('ID')

set.seed(as.numeric(argv[6]))
pdf(file=paste0("./ROC/",modelname,"_ROC.pdf"), width=5, height=5) 
  for(i in 1:length(rs)){
    rocdata=rs[[i]]
    roc1=roc(rocdata$group, as.numeric(rocdata$RS))
    ci1=ci.auc(roc1, method="bootstrap")
    ciVec=as.numeric(ci1)
    if(i == 1){
    	plot(roc1, print.auc=F, col=i+1, legacy.axes=T, main=modelname)
	text(0.3, 0.2-0.05*i, cex=0.6,paste("Train Set",paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }else{
        plot(roc1, print.auc=F, col=i+1, legacy.axes=T, add=T)
    	text(0.3, 0.2-0.05*i, cex=0.6,paste(paste("Test",i-1,sep=""),paste(paste("AUC=", round(ciVec[2]*100, 2), "%"), paste("95% CI:", round(ciVec[1]*100, 2), "%-", round(ciVec[3]*100, 2), "%"),sep="; "),sep=" : "), col=i+1)
    }
  }
  dev.off()

  AUCs$Model <- modelname
  result <- rbind(result,AUCs)

write.table(result,paste0("./ROC/",modelname,"_AUC_result"),sep="\t",quote=F,col.names=NA)



