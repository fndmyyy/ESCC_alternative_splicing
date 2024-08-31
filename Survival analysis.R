library(dplyr)
library(readxl)
library(tidyr)
library(tidyverse)
library(survival)
library(survminer)
library(cutoff)
library(readxl)
rm(list=ls()) 
#data <-read.csv("免疫组化生存分析数据.csv",check.names = F)
data<-as.data.frame(read_xlsx("20240807免疫组化(1).xlsx",sheet = "survival"))

res.cut <- cutoff::logrank(data = data,
                           time = "OS.Time",
                           y = "OS", 
                           x = "SF3B2 positive cell rate",
                           cut.numb = 1, #截点个数
                           n.per = 0.01, #分组后每组样本量占总样本量的最小比例
                           y.per = 0.01, #分组后每组中较少结果的最小比例
                           p.cut = 1,
                           round = 5) #保留几位小数
res.cut[order(res.cut$pvalue,decreasing  = F),]

pvalue<-c()



sur.cut <-surv_cutpoint(data,time= 'OS.Time',event = 'OS',variables = "SF3B2 positive cell rate")
sur.cat <- surv_categorize(sur.cut)
sur.cat$Group<-ifelse(as.numeric(sur.cat$`SF3B2 positive cell rate`)>median(as.numeric(sur.cat$`SF3B2 positive cell rate`)),"High","Low")
fit <- survfit(Surv(OS.Time, OS)~sur.cat$Group, data =sur.cat)
ppvalue<-surv_pvalue(fit)
ppvalue["variable"]<-colnames(data[2])
pvalue<-rbind(pvalue,ppvalue)
ggpar(ggsurvplot(fit,data = sur.cat,pval = TRUE),palette=c("red","blue"),legend.title =paste0("high=",fit$n,"  low=",nrow(data)-fit$n),)
ggsave(paste0(fit,"生存曲线.pdf"))

write.csv(pvalue,"p值表格.csv")




surv_pvalue(fit)
