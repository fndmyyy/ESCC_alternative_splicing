######按照基因进行分组药敏对照
for (i in 1:198)
{
  Cam <- as.data.frame(res[,i])
  colnames(Cam) <- "senstivity"
  Cam$Risk <- group
  boxplot=ggboxplot(Cam, x="Risk", y="senstivity", fill="Risk",
                    xlab="Risk",
                    ylab=paste0(colnames(res)[i], " senstivity (IC50)"),
                    legend.title="Risk",
                    palette=c("green", "red")
  )+
    stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("low risk","high risk")),method="t.test")+#???Ӽ???
    pdf(file=paste0(colnames(res)[i], ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}


####1.联合使用3个数据库分析 2.使用49+179+RNActDrug 联合验证 3.预测剪接事件对药物的响应


library(oncoPredict)
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(dplyr)
library(data.table)
library("parallel")

filename<-list.files(pattern = "vst")

for (i in filename)
{
testExpr<-read.csv(i,row.names = 1)
#testExpr<- testExpr[,grep("T",colnames(testExpr))]
testExpr[,1:ncol(data)] <- as.numeric(unlist(testExpr[,c(1:ncol(data))]))
testExpr<-as.matrix(testExpr)

th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
GDSC2_Expr = readRDS(file=file.path("./",'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path("./","GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

GDSC1_Expr = readRDS(file=file.path("./","GDSC1_Expr (RMA Normalized and Log Transformed).rds"))
GDSC1_Res = readRDS(file = file.path("./","GDSC1_Res.rds"))
GDSC1_Res <- exp(GDSC1_Res) 

CTRP2_Expr = readRDS(file=file.path("./","CTRP2_Expr (TPM, not log transformed).rds"))
CTRP2_Res = readRDS(file = file.path("./","CTRP2_Res.rds"))
CTRP2_Res <- exp(CTRP2_Res) 




unlink("calcPhenotype_Output",recursive = T)

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )
GDSC2 <- fread("calcPhenotype_Output\\DrugPredictions.csv")
GDSC2 <- as.data.frame(GDSC2)
rownames(GDSC2) <- GDSC2$V1
GDSC2 <- GDSC2[,-1]

unlink("calcPhenotype_Output",recursive = T)

calcPhenotype(trainingExprData = GDSC1_Expr,
              trainingPtype = GDSC1_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )
GDSC1 <- fread("calcPhenotype_Output\\DrugPredictions.csv")
GDSC1 <- as.data.frame(GDSC1)
rownames(GDSC1) <- GDSC1$V1
GDSC1 <- GDSC1[,-1]

unlink("calcPhenotype_Output",recursive = T)

calcPhenotype(trainingExprData = CTRP2_Expr,
              trainingPtype = CTRP2_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )
CTRP2 <- fread("calcPhenotype_Output\\DrugPredictions.csv")
CTRP2 <- as.data.frame(CTRP2)
rownames(CTRP2) <- CTRP2$V1
CTRP2 <- CTRP2[,-1]

save(GDSC1,GDSC2,CTRP2,file = paste0(Sys.Date(),i,"药敏原始数据.RData"))

PRISM_Expr<-read.csv("prism-training expression matrix.csv",row.names = 1,check.names = F)
PRISM_Res<-read.csv("prism-training AUC matrix.csv",row.names = 1,check.names = F)
calcPhenotype(trainingExprData = PRISM_Expr,
              trainingPtype = PRISM_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )
PRISM <- fread("calcPhenotype_Output\\DrugPredictions.csv")
PRISM <- as.data.frame(PRISM)
rownames(PRISM) <- PRISM$V1
PRISM <- PRISM[,-1]
save(PRISM,file = paste0(Sys.Date(),i,"药敏原始数据.RData"))
########按照基因表达进行相关性分析
load(list.files(pattern = paste0(i,"药敏原始数据.RData")))
library(foreach)
require(doParallel)
n <- detectCores()
registerDoParallel(n)
Expr<-read.csv(i,row.names = 1)
#Expr<- Expr[,grep("T",colnames(Expr))]
GDSC2<-t(GDSC2)
GDSC1<-t(GDSC1)
CTRP2<-t(CTRP2)
PRISM<-t(PRISM)
#######
a<-1:nrow(Expr)
tasks <- split(a, ceiling(seq_along(a)/n))
results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(GDSC2)) {
      g2 = rownames(GDSC2)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(GDSC2[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(GDSC2[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
GDSC2data_cor <- do.call(rbind, results)

results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(GDSC1)) {
      g2 = rownames(GDSC1)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(GDSC1[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(GDSC1[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
GDSC1data_cor <- do.call(rbind, results)


results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(CTRP2)) {
      g2 = rownames(CTRP2)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(CTRP2[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(CTRP2[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
CTRP2data_cor <- do.call(rbind, results)

results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(PRISM)) {
      g2 = rownames(PRISM)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(PRISM[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(PRISM[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
PRISMdata_cor <- do.call(rbind, results)



#save(list = ls(),file = paste0(Sys.Date(),i,"药敏相关性数据.RData"))

GDSC2data_cor$database<-"GDSC2"
GDSC1data_cor$database<-"GDSC1"
CTRP2data_cor$database<-"CTRP2"
PRISMdata_cor$database<-"PRISM"
data<-rbind(GDSC1data_cor,GDSC2data_cor)
data<-rbind(data,CTRP2data_cor)
data<-rbind(data,PRISMdata_cor)
write.csv(data,paste0(Sys.Date(),i,"药敏相关性数据.RData"))






#data<-fread("10.26 带N样本PRISM 芯片药敏预测，三数据库相关性数据.csv")




gdsc2data<-data.frame(IC50mean=apply(GDSC2,1,mean))
gdsc2data$id<-paste0("GDSC2-",row.names(gdsc2data))

gdsc1data<-data.frame(IC50mean=apply(GDSC1,1,mean))
gdsc1data$id<-paste0("GDSC1-",row.names(gdsc1data))

ctrp2data<-data.frame(IC50mean=apply(CTRP2,1,mean))
ctrp2data$id<-paste0("CTRP2-",row.names(ctrp2data))

prismdata<-data.frame(IC50mean=apply(PRISM,1,mean))
prismdata$id<-paste0("PRISM-",row.names(prismdata))

fourdata<-rbind(gdsc2data,gdsc1data)
fourdata<-rbind(fourdata,ctrp2data)
fourdata<-rbind(fourdata,prismdata)

data$id<-paste0(data$database,"-",data$gene_name2)
data<-merge(data,fourdata,by="id")


write.csv(data,paste0(Sys.Date(),i,"药敏IC50数据.RData"))
}


############################################

write.csv(data,"10.19 芯片表达药敏相关性数据.csv")


#####

sf<-read.csv("剪接因子列表.csv")
sf<-sf[[1]]
sf<-data.frame(RNA.molecule=sf)

data<-fread("RNAactDrug.txt")
data<-merge(data,sf,"RNA.molecule")

############################################# 49对RNA测序数据

library(oncoPredict)
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(dplyr)
library(data.table)


rm(list = ls())
GDSC2_Expr = readRDS(file=file.path("./",'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path("./","GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

GDSC1_Expr = readRDS(file=file.path("./","GDSC1_Expr (RMA Normalized and Log Transformed).rds"))
GDSC1_Res = readRDS(file = file.path("./","GDSC1_Res.rds"))
GDSC1_Res <- exp(GDSC1_Res) 

CTRP2_Expr = readRDS(file=file.path("./","CTRP2_Expr (TPM, not log transformed).rds"))
CTRP2_Res = readRDS(file = file.path("./","CTRP2_Res.rds"))
CTRP2_Res <- exp(CTRP2_Res) 

testExpr<-read.csv("5.4 vst转化去批次效应后的count矩阵.csv",row.names = 1)
#testExpr<- testExpr[,grep("T",colnames(testExpr))]
testExpr<-as.matrix(testExpr)

file.remove("./calcPhenotype_Output/DrugPredictions.csv")
unlink("./calcPhenotype_Output/DrugPredictions.csv", recursive = T, force = T)

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )

cexuGDSC2 <- fread("calcPhenotype_Output\\DrugPredictions.csv")
cexuGDSC2 <- as.data.frame(cexuGDSC2)
rownames(cexuGDSC2) <- cexuGDSC2$V1
cexuGDSC2 <- cexuGDSC2[,-1]

file.remove("./calcPhenotype_Output/DrugPredictions.csv")
unlink("./calcPhenotype_Output/DrugPredictions.csv", recursive = T, force = T)

calcPhenotype(trainingExprData = GDSC1_Expr,
              trainingPtype = GDSC1_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )
cexuGDSC1 <- fread("calcPhenotype_Output\\DrugPredictions.csv")
cexuGDSC1 <- as.data.frame(cexuGDSC1)
rownames(cexuGDSC1) <- cexuGDSC1$V1
cexuGDSC1 <- cexuGDSC1[,-1]

file.remove("./calcPhenotype_Output/DrugPredictions.csv")
unlink("./calcPhenotype_Output/DrugPredictions.csv", recursive = T, force = T)





calcPhenotype(trainingExprData = CTRP2_Expr,
              trainingPtype = CTRP2_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )
cexuCTRP2 <- fread("calcPhenotype_Output\\DrugPredictions.csv")
cexuCTRP2 <- as.data.frame(cexuCTRP2)
rownames(cexuCTRP2) <- cexuCTRP2$V1
cexuCTRP2 <- cexuCTRP2[,-1]
save(cexuCTRP2,file = "10.19 带N样本 49测序 CTRP2（去批次）药敏预测原始数据.RData")
file.remove("./calcPhenotype_Output/DrugPredictions.csv")
unlink("./calcPhenotype_Output/DrugPredictions.csv", recursive = T, force = T)

save(cexuGDSC1,cexuGDSC2,cexuCTRP2,file = "10.19 带N样本 49测序（去批次）药敏预测原始数据.RData")
library(oncoPredict)
library(data.table)
setwd("F:/Document/Oncopredict/10.24 prism")
testExpr<-read.csv("5.4 vst转化去批次效应后的count矩阵.csv",row.names = 1)
testExpr<-as.matrix(testExpr)
PRISM_Expr<-read.csv("prism-training expression matrix.csv",row.names = 1,check.names = F)
PRISM_Res<-read.csv("prism-training AUC matrix.csv",row.names = 1,check.names = F)
PRISM_Expr<-as.matrix(PRISM_Expr)
PRISM_Res<-as.matrix(PRISM_Res)
calcPhenotype(trainingExprData = PRISM_Expr,
              trainingPtype = PRISM_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              printOutput=TRUE,
              minNumSamples = 10, 
              removeLowVaringGenesFrom = 'rawData' )
cexuPRISM <- fread("calcPhenotype_Output\\DrugPredictions.csv")
cexuPRISM <- as.data.frame(cexuPRISM)
rownames(cexuPRISM) <- cexuPRISM$V1
cexuPRISM <- cexuPRISM[,-1]
save(cexuPRISM,file = "10.26 带N样本 49测序 PRISM（去批次）药敏预测原始数据.RData")


######按照基因表达进行相关性分析
load("10.19 带N样本 49测序（去批次）药敏预测原始数据.RData")
load("10.26 带N样本 49测序 PRISM（去批次）药敏预测原始数据.RData")
load("D:/BaiduNetdiskDownload/Document/Oncopredict/10.19 药敏预测/10.19 带N样本 49测序 CTRP2（去批次）药敏预测原始数据.RData")
library(foreach)
require(doParallel)
n <- 48
registerDoParallel(n)
Expr<-read.csv("5.4 vst转化去批次效应后的count矩阵.csv",row.names = 1)
#Expr<- Expr[,grep("T",colnames(Expr))]

GDSC2<-t(cexuGDSC1)
GDSC1<-t(cexuGDSC2)
CTRP2<-t(cexuCTRP2)
PRISM<-t(cexuPRISM)
#######
a<-1:nrow(Expr)
tasks <- split(a, ceiling(seq_along(a)/n))
results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(GDSC2)) {
      g2 = rownames(GDSC2)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(GDSC2[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(GDSC2[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
GDSC2data_cor <- do.call(rbind, results)

results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(GDSC1)) {
      g2 = rownames(GDSC1)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(GDSC1[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(GDSC1[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
GDSC1data_cor <- do.call(rbind, results)


results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(CTRP2)) {
      g2 = rownames(CTRP2)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(CTRP2[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(CTRP2[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
CTRP2data_cor <- do.call(rbind, results)



results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(PRISM)) {
      g2 = rownames(PRISM)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(PRISM[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(PRISM[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
PRISMdata_cor <- do.call(rbind, results)
GDSC2data_cor$database<-"GDSC2"
GDSC1data_cor$database<-"GDSC1"
CTRP2data_cor$database<-"CTRP2"
PRISMdata_cor$database<-"PRISM"
data<-rbind(GDSC1data_cor,GDSC2data_cor)
data<-rbind(data,CTRP2data_cor)
data<-rbind(data,PRISMdata_cor)

save(list = ls(),file = "10.26 带N样本PRISM 49对测序表达 药敏预测(原始和相关性).RData")


write.csv(data,"10.26 带N样本PRISM 49对测序表达 药敏预测，三数据库相关性数据.csv")






data<-fread("10.26 带N样本PRISM 49对测序表达 药敏预测，三数据库相关性数据.csv")
data<-filter(data,pvalue<0.05)



gdsc2data<-data.frame(IC50mean=apply(GDSC2,1,mean))
gdsc2data$id<-paste0("GDSC2-",row.names(gdsc2data))

gdsc1data<-data.frame(IC50mean=apply(GDSC1,1,mean))
gdsc1data$id<-paste0("GDSC1-",row.names(gdsc1data))

ctrp2data<-data.frame(IC50mean=apply(CTRP2,1,mean))
ctrp2data$id<-paste0("CTRP2-",row.names(ctrp2data))

prismdata<-data.frame(IC50mean=apply(PRISM,1,mean))
prismdata$id<-paste0("PRISM-",row.names(prismdata))

fourdata<-rbind(gdsc2data,gdsc1data)
fourdata<-rbind(fourdata,ctrp2data)
fourdata<-rbind(fourdata,prismdata)

data$id<-paste0(data$database,"-",data$gene_name2)
data<-merge(data,fourdata,by="id")

write.csv(data,"10.26 带N样本PRISM 49测序表达药敏相关性数据.csv")

#sf<-read.csv("剪接因子列表.csv")
#sf<-sf[[1]]
#sf<-data.frame(gene_name1=sf)
#sf<-merge(data,sf,by="gene_name1")

#基于剪接事件进行计算

library(foreach)
require(doParallel)
n <- 48
registerDoParallel(n)
Expr<-read.csv("7.25 3699剪接纯原始数据.csv",row.names = 1)
#Expr<- Expr[,grep("T",colnames(Expr))]

GDSC2<-t(cexuGDSC1)
GDSC1<-t(cexuGDSC2)
CTRP2<-t(cexuCTRP2)
PRISM<-t(cexuPRISM)

GDSC2<- GDSC2[,c(grep("T",colnames(GDSC2)),grep("N",colnames(GDSC2)))]
GDSC1<- GDSC1[,c(grep("T",colnames(GDSC1)),grep("N",colnames(GDSC1)))]
CTRP2<- CTRP2[,c(grep("T",colnames(CTRP2)),grep("N",colnames(CTRP2)))]
PRISM<- PRISM[,c(grep("T",colnames(PRISM)),grep("N",colnames(PRISM)))]





#######
a<-1:nrow(Expr)
tasks <- split(a, ceiling(seq_along(a)/n))
results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(GDSC2)) {
      g2 = rownames(GDSC2)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(GDSC2[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(GDSC2[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
GDSC2data_cor <- do.call(rbind, results)

results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(GDSC1)) {
      g2 = rownames(GDSC1)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(GDSC1[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(GDSC1[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
GDSC1data_cor <- do.call(rbind, results)


results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(CTRP2)) {
      g2 = rownames(CTRP2)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(CTRP2[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(CTRP2[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
CTRP2data_cor <- do.call(rbind, results)



results <- foreach(m = tasks) %dopar% {
  
  # 在每个线程维护自己的对象
  gene_name1 <- c() 
  gene_name2 <- c()
  cor_r <- c()
  pvalue <- c()
  for (i in m) {
    g1 = rownames(Expr)[i]
    for (j in 1:nrow(PRISM)) {
      g2 = rownames(PRISM)[j]  
      c_r = cor(as.numeric(Expr[i,]),as.numeric(PRISM[j,]),method="spearman",use="pairwise.complete.obs")
      p = cor.test(as.numeric(Expr[i,]),as.numeric(PRISM[j,]),method ="spearman",use="pairwise.complete.obs")[[3]]
      gene_name1 <- c(gene_name1, g1)
      gene_name2 <- c(gene_name2, g2)
      cor_r <- c(cor_r, c_r)
      pvalue <- c(pvalue, p)
    }
  }
  
  data.frame(gene_name1, gene_name2, cor_r, pvalue)  
}
PRISMdata_cor <- do.call(rbind, results)
GDSC2data_cor$database<-"GDSC2"
GDSC1data_cor$database<-"GDSC1"
CTRP2data_cor$database<-"CTRP2"
PRISMdata_cor$database<-"PRISM"
data<-rbind(GDSC1data_cor,GDSC2data_cor)
data<-rbind(data,CTRP2data_cor)
data<-rbind(data,PRISMdata_cor)

save(list = ls(),file = "10.26 带N样本PRISM 49对测序剪接 药敏预测(原始和相关性).RData")

write.csv(data,"10.26 带N样本PRISM 49对测序剪接药敏预测，三数据库相关性数据.csv")





library(data.table)
library(dplyr)
data<-fread("10.26 带N样本PRISM 49对测序剪接药敏预测，三数据库相关性数据.csv")
data<-filter(data,pvalue<0.05)



gdsc2data<-data.frame(IC50mean=apply(GDSC2,1,mean))
gdsc2data$id<-paste0("GDSC2-",row.names(gdsc2data))

gdsc1data<-data.frame(IC50mean=apply(GDSC1,1,mean))
gdsc1data$id<-paste0("GDSC1-",row.names(gdsc1data))

ctrp2data<-data.frame(IC50mean=apply(CTRP2,1,mean))
ctrp2data$id<-paste0("CTRP2-",row.names(ctrp2data))

prismdata<-data.frame(IC50mean=apply(PRISM,1,mean))
prismdata$id<-paste0("PRISM-",row.names(prismdata))

fourdata<-rbind(gdsc2data,gdsc1data)
fourdata<-rbind(fourdata,ctrp2data)
fourdata<-rbind(fourdata,prismdata)

data$id<-paste0(data$database,"-",data$gene_name2)
data<-merge(data,fourdata,by="id")

write.csv(data,"10.26 带N样本PRISM 49测序剪接药敏相关性数据.csv")




#### 筛选
library(data.table)
library(dplyr)
setDTthreads(48)
cexujianjieyuce<-fread("10.19 49测序剪接药敏相关性数据.csv",check.names = T)
cexubiaodayuce<-fread("10.19 49测序表达药敏相关性数据.csv",check.names = T)
xinpianbiaodayuce<-fread("10.19 芯片表达药敏相关性数据.csv",check.names = T)
cexujianjieyuce<-cexujianjieyuce[which(cexujianjieyuce$IC50meanT<10),]
cexubiaodayuce<-cexubiaodayuce[which(cexubiaodayuce$IC50meanT<10),]
xinpianbiaodayuce<-xinpianbiaodayuce[which(xinpianbiaodayuce$IC50meanT<10),]

cexubiaodayuce<-cexubiaodayuce[,-1]
xinpianbiaodayuce<-xinpianbiaodayuce[,-1]


table(cexubiaodayuce$database)
table(cexujianjieyuce$database)
table(xinpianbiaodayuce$database)
xinpianbiaodayuce<-filter(xinpianbiaodayuce,pvalue<0.05)

sf<-read.csv("剪接因子列表.csv")
sf<-sf[[1]]
sf<-data.frame(gene_name1=sf)
sfcexubiaodayuce<-merge(cexubiaodayuce,sf,by="gene_name1")
sfxinpianbiaodayuce<-merge(xinpianbiaodayuce,sf,by="gene_name1")
table(sfcexubiaodayuce$database)
table(sfxinpianbiaodayuce$database)


write.csv(cexujianjieyuce,"9.7 筛选后 测序剪接药敏相关性数据.csv")
write.csv(sfcexubiaodayuce,"9.7 筛选后 剪接因子测序表达药敏相关性数据.csv")
write.csv(sfxinpianbiaodayuce,"9.7 筛选后 剪接因子芯片表达接药敏相关性数据.csv")


#####t-test and FC




t.test(x1,x2,alternative = "two.sided",paired=T)
FC<-exp(log2(a)-log2(b))
log2FC<-log2(a)-log2(b)




