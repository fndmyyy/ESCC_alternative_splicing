library(ggplot2)
library(survminer)
library(tidyverse)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(limma)
library("ggpubr")
library(tidyr)
library(ggpmisc)
library(readxl)



fileinput

for (filename in fileinput)
{#filename<-""PDL1-SF3B1 mrna.csv"
  method="pearson"
data<-read.csv(filename,check.names = F)

c_r = cor(as.numeric(data[,colnames(data)[1]]),as.numeric(data[,colnames(data)[2]]),method=method,use="pairwise.complete.obs")
p = cor.test(as.numeric(data[,colnames(data)[1]]),as.numeric(data[,colnames(data)[2]]),method =method,use="pairwise.complete.obs")[[3]]

ggplot(data,aes(x=as.numeric(data[,colnames(data)[1]]),y=as.numeric(data[,colnames(data)[2]])))+
  geom_smooth(method="lm",  # 拟合曲线
              se=T,         # 是否添加置信区间
              color="#FA8072")+
  geom_point(colour = "black",size = 2)+             # 散点图
  annotate('text',
    x=min(data[,colnames(data)[1]]),y=max(data[,colnames(data)[2]]),            # 坐标位置，左上角
    label = paste0(method,"'s r=",format(c_r, scientific = FALSE, digits = 2),                                   # 注释标记内容相关性+P值
    "\nP=",format(p, scientific = FALSE, digits = 2)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
  theme_gray()+
  xlab(colnames(data)[1])+
  ylab(colnames(data)[2])
ggsave(paste0(Sys.Date(),filename,"相关性图.pdf"),width = 8,height = 7)}


