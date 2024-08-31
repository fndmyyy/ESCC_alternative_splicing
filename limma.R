setwd("D:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\蛋白质")
library(limma)
library(readxl)
gene<-read.csv("9300蛋白质表达.csv",row.names = 1)
group<-rep(c("N","T"),each=124)
design <- model.matrix(~0+factor(group))
View(design)
colnames(design)=levels(factor(group))
ddd<-gene
rownames(design)=colnames(ddd)
contrast.matrix<-makeContrasts(T-N,levels = design)
contrast.matrix
fit <-lmFit(gene,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
tempOutput = topTable(fit2, coef=1, n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
head(nrDEG)
View(nrDEG)
write.csv(nrDEG,"nrDEG.csv")




####
#nomalizearrarys gene expression.csv
sf<-read.csv("gsesfdeg.csv",header = T)
sf<-sf$genesymbol
geeee<-read.csv("nomalizearrarys gene expression.csv",row.names = 1)
data<-geeee
genesf<-c()
for (i in sf)
{genesf<-rbind(genesf,data[which(row.names(data)==i),])}
write.csv(genesf,"(batch)gsesf表达.csv")

##
data<-read.csv("GSE53625_expression.csv",row.names = 1)
data1<-normalizeBetweenArrays(data)
boxplot(data1)
write.csv(data1,"nomalizearrarys gene expression.csv")
