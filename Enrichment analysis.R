rm(list = ls())
library(readxl)
library(WebGestaltR)
library(data.table)
library(dplyr)
###########剪接因子-剪接事件网络
geneFile0 <-read_xlsx("5.10 蛋白编码WGCNA数据表.xlsx",sheet = "geneInfo")
geneFile0<-geneFile0[,c(2,3)]
refFile <-as.vector(read_xlsx("参考基因列表.xlsx",sheet = 1))

enrich<-c()
# 循环处理每个模块颜色
for (i in as.character(unique(geneFile0$moduleColor))) {
  geneFile <- filter(geneFile0, moduleColor == i)
  geneFile <- geneFile[1]
  colnames(geneFile) <- "ID"
  geneFile <- as.vector(geneFile)
  
  # 尝试运行WebGestaltR函数，捕获错误并跳过该次循环
  tryCatch({
    enrichResult_HALLMARK <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
                                         enrichDatabase = "HALLMARK", enrichDatabaseFile = "h.all.v2022.1.Hs.symbols.gmt", enrichDatabaseType = "genesymbol",
                                         interestGene = geneFile[["ID"]], interestGeneType = "genesymbol",
                                         referenceGene = refFile[["ID"]], referenceGeneType = "genesymbol",
                                         sigMethod = "top", topThr = 50, isOutput = TRUE, outputDirectory = getwd(),
                                         projectName = "HALLMARK", nthreas = 16)
    
    enrichResult_HALLMARK$description <- enrichResult_HALLMARK$geneSet
    enrichResult_HALLMARK$geneSet <- "HALLMARK"
    enrichResult_HALLMARK$moduleColor <- i
    enrich <- rbind(enrich, enrichResult_HALLMARK)
    
  }, error = function(e) {
    # 记录错误
    print(paste("Error in loop", i))  
    
    # 设置标志位,在for循环外部判断是否跳过
    skipLoop <<- TRUE 
  })
  
  # 在for循环外检查标志位 
  if(exists("skipLoop") && skipLoop){
    skipLoop <<- FALSE  
    next
  }
  
}




enrichResult_GOBP <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                 enrichDatabase="geneontology_Biological_Process", interestGene=geneFile[["ID"]], 
                                 interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                 referenceGeneType="genesymbol",sigMethod="top",isOutput=TRUE, outputDirectory=getwd(), 
                                 projectName="GO_BP",nthreas=16)
enrichResult_GOBP$geneSet<-"GO_BP"

enrichResult_KEGG <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                 enrichDatabase="pathway_KEGG", interestGene=geneFile[["ID"]], 
                                 interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                 referenceGeneType="genesymbol",sigMethod="top",isOutput=TRUE, outputDirectory=outputDirectory, 
                                 projectName="KEGG",nthreas=16)

enrichResult_KEGG$geneSet<-"KEGG"


enrichResult_HALLMARK <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                 enrichDatabase="HALLMARK",enrichDatabaseFile = "h.all.v2022.1.Hs.symbols.gmt",enrichDatabaseType="genesymbol",interestGene=geneFile[["ID"]], 
                                 interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                 referenceGeneType="genesymbol",sigMethod="top",topThr=10,isOutput=TRUE, outputDirectory=outputDirectory, 
                                 projectName="HALLMARK",nthreas=16)
enrichResult_HALLMARK$description<-enrichResult_HALLMARK$geneSet
enrichResult_HALLMARK$geneSet<-"HALLMARK"


A<-rbind(enrichResult_GOBP,enrichResult_KEGG,enrichResult_HALLMARK)
#绘图
library(ggplot2)
A[which(A$pValue==0),"pValue"]=2.2E-16
A$LogP <- -log10(A$pValue)
x <- ggplot(A,aes(enrichmentRatio,description))+
  geom_bar(aes(y=reorder(description,enrichmentRatio),x=enrichmentRatio,fill=LogP)
           ,stat='identity',width=0.7)+
  scale_fill_gradient(low="#FFCC33",high="#CC6666")+
  theme_bw()+facet_grid(geneSet~.,scales = "free",space = "free")+ 
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  ggtitle(paste0(richgene,"调控富集"))

ggsave(x,filename = paste0(richgene,"调控三数据库富集结果.pdf"),width = 9,height = 10,dpi=300)
write.csv(A,paste0(richgene,"调控三数据库富集结果.csv"))

###########lnc调控
total<-c()
outputDirectory<-getwd()
geneFile1 <-fread("6.23 lncrna剪接相关性+结合预测(lncrna-mrna).csv")
refFile <-as.vector(read_xlsx("参考基因列表.xlsx",sheet = 1))
richgene1<-unique(geneFile1$geneName)
for (richgene in richgene1)
{geneFile<-geneFile1[which(geneFile1$geneName==richgene),]
geneFile<-as.data.frame(geneFile$pairGeneName)
colnames(geneFile)<-"ID"
geneFile<-as.vector(geneFile)

enrichResult_GOBP <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                 enrichDatabase="geneontology_Biological_Process", interestGene=geneFile[["ID"]], 
                                 interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                 referenceGeneType="genesymbol",sigMethod="top",isOutput=TRUE, outputDirectory=outputDirectory, 
                                 projectName="GO_BP",nthreas=16)
enrichResult_GOBP$geneSet<-"GO_BP"


enrichResult_KEGG <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                 enrichDatabase="pathway_KEGG", interestGene=geneFile[["ID"]], 
                                 interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                 referenceGeneType="genesymbol",sigMethod="top",isOutput=TRUE, outputDirectory=outputDirectory, 
                                 projectName="KEGG",nthreas=16)

enrichResult_KEGG$geneSet<-"KEGG"


enrichResult_HALLMARK <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                    enrichDatabase="HALLMARK",enrichDatabaseFile = "h.all.v2022.1.Hs.symbols.gmt",enrichDatabaseType="genesymbol",interestGene=geneFile[["ID"]], 
                                     interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                     referenceGeneType="genesymbol",sigMethod="top",topThr=10,isOutput=TRUE, outputDirectory=outputDirectory, 
                                    projectName="HALLMARK",nthreas=16)
enrichResult_HALLMARK$description<-enrichResult_HALLMARK$geneSet
enrichResult_HALLMARK$geneSet<-"HALLMARK"

A<-rbind(enrichResult_GOBP,enrichResult_KEGG,enrichResult_HALLMARK)
A$incgene<-richgene0
total<-rbind(total,A)
}
#########
for (richgene in richgene1) {
  tryCatch({
    geneFile<-geneFile1[which(geneFile1$geneName==richgene),]
    geneFile<-as.data.frame(geneFile$pairGeneName)
    colnames(geneFile)<-"ID"
    geneFile<-as.vector(geneFile)
    
    enrichResult_GOBP <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                     enrichDatabase="geneontology_Biological_Process", interestGene=geneFile[["ID"]], 
                                     interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                     referenceGeneType="genesymbol",sigMethod="top",isOutput=TRUE, outputDirectory=outputDirectory, 
                                     projectName="GO_BP",nthreas=16)
    enrichResult_GOBP$geneSet<-"GO_BP"
    
    
    enrichResult_KEGG <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                     enrichDatabase="pathway_KEGG", interestGene=geneFile[["ID"]], 
                                     interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                     referenceGeneType="genesymbol",sigMethod="top",isOutput=TRUE, outputDirectory=outputDirectory, 
                                     projectName="KEGG",nthreas=16)
    
    enrichResult_KEGG$geneSet<-"KEGG"
    
    
    enrichResult_HALLMARK <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                         enrichDatabase="HALLMARK",enrichDatabaseFile = "h.all.v2022.1.Hs.symbols.gmt",enrichDatabaseType="genesymbol",interestGene=geneFile[["ID"]], 
                                         interestGeneType="genesymbol", referenceGene=refFile[["ID"]], 
                                         referenceGeneType="genesymbol",sigMethod="top",topThr=10,isOutput=TRUE, outputDirectory=outputDirectory, 
                                         projectName="HALLMARK",nthreas=16)
    enrichResult_HALLMARK$description<-enrichResult_HALLMARK$geneSet
    enrichResult_HALLMARK$geneSet<-"HALLMARK"
    A<-rbind(enrichResult_GOBP,enrichResult_KEGG,enrichResult_HALLMARK) 
    A$incgene<-richgene
    total<-rbind(total,A)
  }, error = function(e) {
    message("Error encountered, skipping gene:", richgene)
  })
}
#绘图
library(ggplot2)
A[which(A$pValue==0),"pValue"]=2.2E-16
A$LogP <- -log10(A$pValue)
x <- ggplot(A,aes(enrichmentRatio,description))+
  geom_bar(aes(y=reorder(description,enrichmentRatio),x=enrichmentRatio,fill=LogP)
           ,stat='identity',width=0.7)+
  scale_fill_gradient(low="#FFCC33",high="#CC6666")+
  theme_bw()+facet_grid(geneSet~.,scales = "free",space = "free")+ 
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  ggtitle(paste0(richgene,"调控富集"))

ggsave(x,filename = paste0(richgene,"调控三数据库富集结果.pdf"),width = 9,height = 10,dpi=300)
write.csv(A,paste0(richgene,"调控三数据库富集结果.csv"))



############
library(readxl)
library(ggplot2)
library(showtext)
inform<-excel_sheets("Fig 1 A-F对应原始富集后数据.xlsx")
for (i in inform)
{data<-read_xlsx("Fig 1 A-F对应原始富集后数据.xlsx",sheet = i)
p <- ggplot(data,aes(x=enrichmentRatio,y=description,colour=-1*log10(pValue),size=-1*log10(pValue)))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "#FFFACD",high = "#01665e")+
  theme_bw()+
  theme(axis.title = element_text(
    family = "sans",##坐标轴标签字体
    face='plain', ##字体外形（粗斜体等） ##字体大小
    lineheight = 1),##标签行间距的倍数
    axis.text = element_text(
      family = "sans",##字体
      face="plain", ##字体外形（粗斜体等）
      color="black"
      ))+ 
  ylab("Pathway")+
  xlab("enrichmentRatio")+
  labs(color=expression(-log[10](pValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="plain",color="black",angle=0,vjust=1))
plot <- p+theme_bw()
ggsave(paste0(i,"富集气泡图.pdf"),height =70 ,width =188,units = "mm" )
}
p




