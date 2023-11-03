#DNA methylation analysis
#DNA甲基化芯片数据分析


#####1.安装依赖包##########

rm(list = ls())   
options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror


# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGG.db",ask = F,update = F)


BiocManager::install("minfi",ask = F,update = F)
BiocManager::install("ChAMP",ask = F,update = F)
BiocManager::install("methylationArrayAnalysis",ask = F,update = F) 
BiocManager::install("wateRmelon",ask = F,update = F) 

BiocManager::install(c("GSEABase","GSVA","clusterProfiler" ),ask = F,update = F)
BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
BiocManager::install(c("org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)

options()$repos
install.packages('WGCNA')
install.packages(c("FactoMineR", "factoextra"))
install.packages(c("ggplot2", "pheatmap","ggpubr"))
library("FactoMineR")
library("factoextra")

library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)


#####2.芯片数据分析##########

#使用GEO数据库下载的甲基化信号值矩阵文件

require(GEOquery)
require(Biobase)

eset <- getGEO("GSE68777",destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])

#获取临床信息
pD.all <- pData(eset[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)
names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)


library(ChAMP)
# beta 信号值矩阵里面不能有NA值
myLoad=champ.filter(beta = beta.m ,pd = pD)
myLoad
save(myLoad,file = 'step1-output.Rdata')


#######3.质量控制##########

rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file <- 'step1-output.Rdata')

#耗时步骤，运行一次后，就注释掉
if(F){
  myNorm <- champ.norm(beta=myLoad$beta,arraytype <-"450K",cores=5)
  dim(myNorm) 
  pD <- myLoad$pd
  save(myNorm,pD,file <- 'step2-champ_myNorm.Rdata')
}


load(file = 'step2-champ_myNorm.Rdata')
# 原来的450K经过质控过滤后是400K啦
beta.m <- myNorm
group_list <- myLoad$pd$group
groupdim(beta.m) 
# 下面是表达矩阵标准3张图质量控制手段，生信技能树原创
if(T){
  dat <- t(beta.m)
  dat[1:4,1:4] 
  library("FactoMineR")#画主成分分析图需要加载这两个包
  library("factoextra")  
  # 因为甲基化芯片是450K或者850K，几十万行的甲基化位点，所以PCA不会太快
  dat.pca <- PCA(dat , graph = FALSE) 
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave('all_samples_PCA.png')
  
  dat=beta.m
  dat[1:4,1:4] 
  cg=names(tail(sort(apply(dat,1,sd)),1000))#apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
  n=t(scale(t(dat[cg,]))) # 'scale'可以对log-ratio数值进行归一化
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(n)  
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,filename = 'heatmap_top1000_sd.png')
  dev.off()
  
  exprSet=beta.m
  pheatmap::pheatmap(cor(exprSet)) 
  # 组内的样本的相似性应该是要高于组间的！
  colD=data.frame(group_list=group_list)
  rownames(colD)=colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet),
                     annotation_col = colD,
                     show_rownames = F,
                     filename = 'cor_all.png')
  dev.off() 
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  # M=cor(log2(exprSet+1)) 
  M=cor(exprSet)
  pheatmap::pheatmap(M,annotation_col = colD)
  pheatmap::pheatmap(M,
                     show_rownames = F,
                     annotation_col = colD,
                     filename = 'cor_top500.png')
  dev.off() 
  
}


myLoad 
beta.m=myLoad$beta
# 下面是使用 wateRmelon 进行 归一化代码，目前被主流抛弃
if(F){
  # 使用 wateRmelon 进行 归一化
  library("wateRmelon")
  beta.m=beta.m[rowMeans(beta.m)>0.005,]
  pdf(file="rawBox.pdf")
  boxplot(beta.m,col = "blue",xaxt = "n",outline = F)
  dev.off()
  beta.m = betaqn(beta.m)
  pdf(file="normalBox.pdf")
  boxplot(beta.m,col = "red",xaxt = "n",outline = F)
  dev.off()
  
  # 然后进行简单的QC
  group_list=myLoad$pd$group
  pdf(file="densityBeanPlot.pdf")
  par(oma=c(2,10,2,2))
  densityBeanPlot(beta.m, sampGroups = group_list)
  dev.off()
  pdf(file="mdsPlot.pdf")
  mdsPlot(beta.m, numPositions = 1000, sampGroups = group_list)
  dev.off()
  
  # 后续针对 beta.m 进行差异分析, 比如 minfi 包
  grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
  M = getM(grset)
  # 因为甲基化芯片是450K或者850K，几十万行的甲基化位点，统计检验通常很慢。
  dmp <- dmpFinder(M, pheno=group_list, type="categorical")
  dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
  dim(dmpDiff)
}







########DEM-champ##########


library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')
myLoad    # 存储了甲基化信号矩阵和表型信息。
load(file = 'step2-champ_myNorm.Rdata')
group_list=myLoad$pd$group
table(group_list)
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
head(myDMP[[1]])
save(myDMP,file = 'step3-output-myDMP.Rdata')
# 还可以调动交互式界面修改阈值，调整差异化探针
# DMP.GUI(DMP=myDMP[[1]],beta=myNorm,group_list)

# 下面的分析, 非常的消耗计算资源
# 如果你有时间，就折腾
if(F){
  myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
  DMR.GUI(DMR=myDMR)
  
  myBlock <- champ.Block(beta = myNorm,pheno=group_list,arraytype="450K")
  head(myBlock$Block)
  Block.GUI(Block=myBlock,beta = myNorm,pheno=group_list,
            runDMP=TRUE,compare.group=NULL,arraytype="450K")
  
  myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]],
                       DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
  
  head(myGSEA$DMP)
  head(myGSEA$DMR)
  myEpiMod <- champ.EpiMod(beta=myNorm,pheno=group_list)
  
}


#########compare-champ-minfi############

library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')
beta.m=myLoad$beta
# 后续针对 beta.m 进行差异分析, 比如 minfi 包
grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
M = getM(grset)
group_list=myLoad$pd$group
# 因为甲基化芯片是450K或者850K，几十万行的甲基化位点，统计检验通常很慢。
dmp <- dmpFinder(M, pheno=group_list, type="categorical")
dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
dim(dmpDiff)

load(file = 'step3-output-myDMP.Rdata')
champDiff=myDMP[[1]]

dim(dmpDiff)
dim(champDiff)
length(intersect(rownames(dmpDiff),rownames(champDiff)))

source('functions_DEM.R')
visual_champ_DEM(myLoad,myDMP,group_list,pro='test')


###########go-kegg##############


library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)

load(file = 'step3-output-myDMP.Rdata')
deg=myDMP[[1]]
head(deg)
length(unique(deg$gene)) 
deg$g=ifelse(abs(deg$logFC) < 0.2,'stable',
             ifelse(deg$logFC > 0.2,'UP','DOWN'))
table(deg$g)
head(deg)
deg$symbol=deg$gene
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] ) 
source('kegg_and_go_up_and_down.R')
# 有时候kegg数据库抽风，这个函数会失败，这个锅，Y叔背
run_kegg(gene_up,gene_down,pro='test_methy')
# 需要多go数据库的3个条目进行3次富集分析，非常耗时。
# 所以我注释掉了下面的代码
# run_go(gene_up,gene_down,pro='npc_VS_normal')

# 下面的GO肯定会慢，但是会成功运行
go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_up_GO_all_barplot.png') 
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_down_GO_all_barplot.png')



##########visualization###############

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

load(file = 'anno_DEG.Rdata') 

gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] ) 
gene_down
gene_up
# 下面的函数依赖于kegg数据库连接网速，有时候中国大陆的朋友会失败
enrichKK <- enrichKEGG(gene         =  gene_up,
                       organism     = 'hsa',
                       #universe     = gene_all,
                       pvalueCutoff = 0.1,
                       qvalueCutoff =0.1)
head(enrichKK)[,1:6] 
browseKEGG(enrichKK, 'hsa04512')
dotplot(enrichKK)
enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
enrichKK 

#(3)可视化
#条带图
par(mfrow=c(2,1))
barplot(enrichKK,showCategory=20)
#气泡图
dotplot(enrichKK)
#下面的图需要映射颜色，设置和示例数据一样的geneList
geneList = deg$logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)
#(3)展示top5通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(enrichKK, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(enrichKK, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#Enrichment Map
emapplot(enrichKK)
#(4)展示通路关系,仅仅是针对于GO数据库结果。
# goplot(enrichKK)
#(5)Heatmap-like functional classification
heatplot(enrichKK,foldChange = geneList)

















