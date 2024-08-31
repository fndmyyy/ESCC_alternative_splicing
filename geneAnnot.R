rm(list = ls())
library(readxl)
library(xlsx)
setwd("E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\5.3基因类型识别")


###
###
allgene<-read_xlsx("E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\5.3基因类型识别\\5.4 基因分类文件.xlsx",sheet = 4)
allgene<-allgene[,c(7,8,9)]
coding<-read_xlsx("E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\5.3基因类型识别\\5.4 基因分类文件.xlsx",sheet = "proteincoding")
coding<-coding[,c(7,8,9)]
noncoding<-read_xlsx("E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\5.3基因类型识别\\5.4 基因分类文件.xlsx",sheet = "noncoding")
noncoding<-noncoding[,c(7,8,9)]
###
genee<-read_xls("410 T 差异表达数据DESeq差异表达数据042223.xls")
dd<-gsub(c("ID"),c("gene_name"),colnames(genee))
colnames(genee)<-dd
####
genetabel<-read.csv("E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\4.22 不删ENSG\\49对表达\\5.4 vst转化去批次效应后的count矩阵.csv")
genedeg<-read.table("E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\4.22 不删ENSG\\49对表达\\5.4 DESeq差异表达数据.txt",check.names = F,header = T)
asevent<-read_xlsx("E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\4.22 不删ENSG\\5.4 49对剪接纯处理数据\\5.4 显著剪接事件列表.xlsx",sheet = 1)
###
genetabel<-merge(genetabel,allgene,by = "gene_name")
genedeg<-merge(genedeg,allgene,by = "gene_name")
genee<-merge(genee,allgene,by = "gene_name")
###
asevent<-merge(asevent,allgene,by = "gene_name")

write.csv(genee,"(5.9注释)410 T 差异表达数据DESeq差异表达数据042223.csv")
write.csv(genetabel,"E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\4.22 不删ENSG\\49对表达\\5.4 (注释基因类型)vst转化去批次效应后的count矩阵.csv")
write.csv(genedeg,"E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\4.22 不删ENSG\\49对表达\\5.4 (注释基因类型)DESeq差异表达数据.csv" )
write.csv(asevent,"E:\\BaiduSyncdisk\\大创下游--全部重做\\(4.4)大创重做\\4.22 不删ENSG\\5.4 49对剪接纯处理数据\\5.4 (注释基因类型)显著剪接事件列表.csv")