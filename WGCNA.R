rm(list = ls()) ## 魔幻操作，一键清空~
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型

###
library(readxl)
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table) #多核读取文件
### 启用WGCNA多核计算
femData = read.csv( paste0(i,"处理后数据.csv") ,header = T,row.names = 1) #载入基因表达量数据
#读取表型数据
allTraits=read.csv(paste0(i,"处理后队列数据.csv"),row.names = 1)
aa<-model.matrix(~0+factor(allTraits[[1]]))
row.names(aa)<-row.names(allTraits)
allTraits<-as.data.frame(aa)
#行为基因，列为不同样本的基因表达量或其他信息
datExpr0 = as.data.frame(t(femData))  #提取加转置
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#表型和样本匹配
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, row.names(allTraits));
datTraits = allTraits[traitRows, ];
rownames(datTraits) = row.names(allTraits);
collectGarbage()
#可视化表型数据与基因表达量数据的联系，重构样本聚类树
sampleTree = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
#用颜色代表关联度
pdf(file="3_Sample_dendrogram_and_trait_heatmap.pdf",width=8 ,height= 6)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
abline(h = 250, col = "red")
dev.off()
clust = cutreeStatic(sampleTree, cutHeight = 225, minSize = 10)
keepSamples = (clust==1)
table(keepSamples)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#表型和样本匹配
#构建表达网络
#构建自动化网络和检测模块(基因，un=12,)
#选择软阈值
allowWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
type = "unsigned"
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="unsigned")
par(mfrow = c(1,2))
cex1 = 0.9
#无标度拓扑拟合指数
pdf(file="3_Sample_dendrogram_and_trait_heatmap.pdf",width=8 ,height= 6)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.8,col="red")  #查看位于0.9以上的点，可以改变高度值
dev.off()
#平均连接度
pdf(file="3_Sample_dendrogram_and_trait_heatmap.pdf",width=8 ,height= 6)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
power=sft$powerEstimate
dev.off()
#连通度验证
pdf(file="3_Sample_dendrogram_and_trait_heatmap.pdf",width=8 ,height= 6)
k <- softConnectivity(datE=datExpr,power=power) 
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()
#net = blockwiseModules(datExpr, power = 5, TOMType = "unsigned", loadTOM = TRUE,saveTOMFileBase = "femaleMouseTOM")
#构建网络，模块检测
allowWGCNAThreads()
enableWGCNAThreads(nThreads = 1*parallel::detectCores()) 

#TOM = TOMsimilarityFromExpr(datExpr, power =5);
TOM = TOMsimilarityFromExpr(datExpr, power = power)
net = blockwiseModules(datExpr, power = power, maxBlockSize = nGenes,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3,loadTOM = TRUE,saveTOMFileBase = "femaleMouseTOM")
save(list=ls(),file=paste0(i,"WGCNA.RData"))

#查看划分的模块数和每个模块里面包含的基因个数
table(net$colors)

#层次聚类树状图
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#保存分配模块和模块包含的基因信息。
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
#模块-表型数据关联
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# 用热图的形式展示相关系数
par(mar=c(3,8,2,1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               cex.lab = 0.5,
               main = paste("Module-trait relationships"))
#colors = greenWhiteRed(50)不适用于红绿色盲患者，建议用 blueWhiteRed代替.
#基因与表型数据的关系、重要模块：基因显著性和模块成员模块与表型数据关联并识别重要基因
tumor= as.data.frame(datTraits$tumor);
names(tumor) = "tumor";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, tumor, use = "p"));#和癌症性状的关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(tumor), sep="");
names(GSPvalue) = paste("p.GS.", names(tumor), sep="");

#模块内分析：鉴定具有高GS和高MM的基因
#使用GS和MM测量，可以识别与体重高度相关的基因，以及感兴趣的模块中的高度相关的成员。这个例子中，体重与棕色模块的关联度较高，因此我们在棕色模块中绘制基因显著性和模块成员关系的散点图。注：
#GS：所有基因表达谱与这个模块的eigengene的相关性（cor）。每一个值代表这个基因与模块之间的关系。如果这个值的绝对值接近0，那么这个基因就不是这个模块中的一部分，如果这个值的绝对值接近1，那么这个基因就与这个模块高度相关。
#MM：基因和表型性状比如体重之间的相关性的绝对值。为了将表型特征信息与共表达网络联合起来，比如体重与哪个模块高度相关。每一个基因的表达值与表型性状之间的相关性的绝对值。0表示这个基因与这个性状不相关，1表示高度相关。如果一个模块中的基因都有这个性状高度相关，那么这个模块也就与这个性状高度相关。
#运行以下代码可视化GS和MM
#基因画了以下四个："turquoise","yellow","skyblue","grey60","green"
module = "grey60"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,abline = T)
#图中的每一个点代表一个基因，横坐标值表示基因与模块的相关性，纵坐标值表示基因与表型性状的相关性，这里可以看出与性状高度显著相关的基因往往是与这个性状显著相关的模块中的重要元素。
#输出网络分析结果，我们在找到了与我们感兴趣的表征高度相关的模块，并通过MM测量确定了核心的参与者（核心基因hubgene），我们现在需要将这些数据合并起来，并输出为结果文件。
names(datExpr)#会返回所有在分析中的基因ID
#返回属于棕色模块的基因ID
names(datExpr)[moduleColors=="brown"]

probes = names(datExpr) # 匹配信息
probes2annot = probes;
geneInfo0 = data.frame(
  geneSymbol = probes,
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue);
#按照与肿瘤的显著水平将模块进行排序:
modOrder = order(-abs(cor(MEs, tumor, use = "p")));
#添加模块成员的信息：
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.tumor));  # 排序
geneInfo = geneInfo0[geneOrder, ]
#输出为CSV格式，可用fix(geneInfo)在R中查看：
write.csv(geneInfo, file = "geneInfo.csv")











#网络可视化
#5.0 参数设置
library(WGCNA)
options(stringsAsFactors = FALSE);
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#5.1 可视化基因网络
#可视化加权网络的方法之一是制作热图。热图的每行每列代表一个基因，浅色代表低邻接（重叠）；深色代表高邻接(重叠)。以下代码是将一步法和逐步法的基础上绘制的热图，不适用于逐块分析法，如需要展示逐步法，需要修改代码将每块block进行可视化。
#计算TOM矩阵
#TOM = TOMsimilarityFromExpr(datExpr, power =12);
dissTOM = 1-TOM;
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
sizeGrWindow(9,9)
library(gplots) # 需要先安装这个包才能加载。
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col=myheatcol)

TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#生成的热图可能需要大量的时间。可以限制基因的数量来加快绘图。但是基因子集的树状图看起来与所有基因的树状图不同，下面随机选取400个基因进行绘图：
nSelect = 5000
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
#改变热图的深色背景为白色背景：
library(gplots) # 需要先安装这个包才能加载。
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)

#5.2 可视化表征基因网络
#研究找到的模块之间的关系，可以使用eigengene表征基因作为代表轮廓，通过特征基因相关性来量化模块的相似性。该包包含一个函数plotEigengeneNetworks，可以生成eigengene网络的摘要图。
# 重新计算模块的eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# 提取肿瘤的表型数据
#tumor = as.data.frame(datTraits$tumor);
tumor= as.data.frame(datTraits$tumor)
paratumor = as.data.frame(datTraits$paratumor )
#names(tumor) = "tumor"
names(tumor) = "tumor"
names(paratumor) = "paratumor"
# 加入到相应的模块
MET = orderMEs(MEs,greyLast = T)
MET = orderMEs(cbind(MEs,tumor,paratumor))

#画图
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(4,4.5,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
#以上结果会生成特征模块与体重数据的聚类图和热图。要想拆分聚类图和热图，可以用以下代码实现。
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(6,6,6,6),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


#从图中结果可知，体重与模块MEbrown、MEred、MEblue的关系更加密切。
#6. 将网络导出到网络可视化软件
#第六步是我们最想要的结果，也是每篇文献中最主要的一个图，就是hub基因的互作关系网络图。这步会告诉你如何将必要的数据导出，以供其他软件进行绘图，例如VisANT、Cytoscape。
#6.0 参数设置与数据导入
getwd();
workingDir = ".";
setwd(workingDir);
library(WGCNA)
options(stringsAsFactors = FALSE);
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
lnames
#6.1 输出到VisANT软件所需的数据
#TOM = TOMsimilarityFromExpr(datExpr, power = 12);
annot = read.csv(file = "GeneAnnotation.csv");
module = "brown";
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )
#因为棕色模块相当大，我们可以严格控制输出的hubgene的个数为30个以内在这个模块中。
nTop = 30;nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )
#以上导出的数据可以用VisANT进行编辑，绘制互作网络。

#6.2 输出到Cytoscape
#Cytoscape 允许用户输入边缘文件和节点文件，允许用户指定例如连接权重和节点颜色。在这里，我们向 Cytoscape 展示了两个模块（红色和棕色模块）的输出。
#TOM = TOMsimilarityFromExpr(datExpr, power =12);
# 选择棕色和红色的模块
#蛋白画了以下四个："lightyellow","midnightblue","tan","turquoise"
#基因画了以下四个："turquoise","yellow","skyblue","grey60"
for (i in c("turquoise","yellow","skyblue","grey60"))
{
  modules = i
  probes = names(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  modGenes = probes[match(modProbes, probes)];
  # 选择相关的 TOM矩阵
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = T,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]);
}







