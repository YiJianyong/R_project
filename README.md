![image](https://github.com/YiJianyong/R_project/assets/78517435/40360fa1-8a9a-4457-bd4b-abf075bbddb1)<!--
 * @Author: your name
 * @Date: 2022-04-16 22:44:40
 * @LastEditTime: 2022-04-16 22:48:48
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: \Markdowne:\Windows-SSD\Program Files (x86)\Common Files\Designer\R_teamwork\R_project\R_project\README.md
-->
# DNA甲基化

## 甲基化芯片数据上游处理

[甲基化芯片的一般分析流程](https://mp.weixin.qq.com/s/JHrL_DqgQY6Yh18vHySKYg)

### 1.数据下载
[数据库中数据规律及处理函数](https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&amp;mid=2247486063&amp;idx=1&amp;sn=156bee5397e979722b36b78284188538&amp;scene=21#wechat_redirect)
+ 下载最原始的芯片数据，根据不同芯片在处理成表达矩阵
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
  + 


+ 质量控制

+ 归一化

## 甲基化芯片数据下游处理  
### 1.差异甲基化
+ 差异化CpG位点：找出一个一个的差异甲基化位点 
+ 差异化CpG区域：找出一个连续不断比较长的差异片段
+ 更大范围的差异化region区域


### 2.可视化分析





### 3.功能富集分析



### 4.相关性分析










