
# DNA甲基化分析

## 甲基化芯片数据上游处理

[甲基化芯片的一般分析流程](https://mp.weixin.qq.com/s/JHrL_DqgQY6Yh18vHySKYg)

[DNA甲基化数据分析专题总结](https://blog.csdn.net/weixin_43569478/article/details/108079702)

### 1.数据下载
[数据库中数据规律及处理函数](https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&amp;mid=2247486063&amp;idx=1&amp;sn=156bee5397e979722b36b78284188538&amp;scene=21#wechat_redirect)
+ 下载最原始的芯片数据，根据不同芯片再处理成表达矩阵
  + 数据库获取：GEO   TCGA
  + 使用GDS号下载的数据需要用**Table()函数**处理得到表达矩阵、用**Meta()函数**得到描述信息、GDS2eSet()函数把它转换成expression set对象
  + 使用GSE号返回的对象就是expression set对象，常用处理函数：geneNames/sampleNames/pData/**exprs**
  + 使用GPL号返回的对象和GDS号一样需要用Table和Meta函数进行处理

```
downGSE <- function(studyID = "GSE1009", destdir = ".") {

    library(GEOquery)
    eSet <- getGEO(studyID, destdir = destdir, getGPL = F)

    exprSet = exprs(eSet[[1]])
    pdata = pData(eSet[[1]])

    write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
    write.csv(pdata, paste0(studyID, "_metadata.csv"))
    return(eSet)

}
```
+ 下载一些文章中作者已经做完归一化后的表达矩阵数据，那就可以少一些步骤


### 2.数据预处理
+ 对原始芯片数据进行处理取决于芯片的平台：这里主要介绍主流的illumina相关芯片
  + 用lumi包处理bead系列表达芯片

```
library(lumi)
fileName <- 'GSE30669_HEK_Sample_Probe_Profile.txt' # Not Run
x.lumi <- lumiR.batch(fileName) ##, sampleInfoFile='sampleInfo.txt')
pData(phenoData(x.lumi))
## Do all the default preprocessing in one step
lumi.N.Q <- lumiExpresso(x.lumi)
### retrieve normalized data
dataMatrix <- exprs(lumi.N.Q)
```
  + 从GEO数据库中直接下载的数据

```
library(GEOquery)
library(limma)
GSE30669 <- getGEO('GSE30669', destdir=".",getGPL = F)
exprSet=exprs(GSE30669[[1]])
GSE30669[[1]]
pdata=pData(GSE30669[[1]])
exprSet=exprs(GSE30669[[1]])
```

+ 对于探针的甲基化水平，常见的定量方式包括beta值和M值两种
  + beta值：**M/(M + U + offset)**  U是代表非甲基化信号强度，M代表甲基化信号强度，offset是偏移量，防止出现分母等于零的现象
  + M值：**log2(M/U)**
  + beta值主要用于差异分析，M值适用于样本间的特征比较

+ p值：p值越小，探针的质量越高，可信度越高

+ 质量控制：主要是探针水平的过滤，包括三个方面
  + 1.**过滤掉可信度较低的探针**：通常认为p值大于0.01的探针信号是不可靠的
  + 2.过滤掉覆盖了SNP位点的探针：探针覆盖了SNP位点可能会影响其杂交情况
  + 3.过滤掉位于性染色体上的探针：样本中混合了两种性别，无法衡量性染色体上的甲基化水平差异

+ 标准化：视情况选择合适的方法
  + 对于用于多样本的差异比较的话：可以用preprocessFunnorm()函数
  + 对于不需要进行差异比较的单样本数据：可以用preprocessQuantile()函数

## 甲基化芯片数据下游处理  
### 1.差异甲基化
+ 差异化CpG位点：找出一个一个的差异甲基化位点 
+ 差异化CpG区域：找出一个连续不断比较长的差异片段
+ 更大范围的差异化region区域


### 2.可视化分析





### 3.功能富集分析



### 4.相关性分析










