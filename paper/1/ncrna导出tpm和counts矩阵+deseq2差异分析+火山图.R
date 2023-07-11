jname1 <- "Control_ncRNA.csv"
jname2 <- "AS_ncRNA.csv"


#设置读取文件路径
path0 <- "D:/data_anlysis/asncRNA"
pathc <- "D:/data_anlysis/asncRNA/Control_Bulk/ncRNA"
pathas <- "D:/data_anlysis/asncRNA/AS_Bulk/ncRNA"


#设置输出文件名称
cname <- "ncrna_control.csv"
cnamesum <- "ncrna_control_sum.csv"
asname <- "ncrna_as.csv"
asnamesum <- "ncrna_as_sum.csv"

cname_count <- "ncrna_control_count.csv"
cnamesum_count <- "ncrna_control_sum_count.csv"

asname_count <- "ncrna_as.csv_count"
asnamesum_count <- "ncrna_as_sum_count.csv"

cinfo <- "cinformation.csv"
asinfo <- "asinformation.csv"

infotestname <- "linkna_allinfo_deseq2.csv"


#使用tximport导出样本tpm表达矩阵
setwd(path0)
library(DESeq2)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)

sample<- dir(file.path(pathc))
files<-file.path(pathc,sample,"abundance.h5")
names(files)<-sample

library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

tx2gene <- t2g  #从biomart包中拉去整体的对应表格（虽然似乎不太全面）

txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)
txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")


write.csv(txi$abundance,cname,quote = F,row.names = T,col.names = T)
write.csv(txi.sum$abundance,cnamesum,quote = F,row.names = T,col.names = T)
write.csv(txi$counts,cname_count,quote = F,row.names = T,col.names = T)
write.csv(txi.sum$counts,cnamesum_count,quote = F,row.names = T,col.names = T)


##########从大矩阵当中抓取含有特定名字的小分类矩阵


#对AS组进行重复的矩阵导出操作
#使用tximport导出样本tpm表达矩阵
setwd(path0)

library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)

sample<- dir(file.path(pathas))
files<-file.path(pathas,sample,"abundance.h5")
names(files)<-sample

library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

tx2gene <- t2g  #从biomart包中拉去整体的对应表格（虽然似乎不太全面）

txi<-tximport(files,type = "kallisto",tx2gene = tx2gene,ignoreAfterBar = TRUE,ignoreTxVersion = TRUE)
txi.kallisto<-tximport(files,type="kallisto",txOut=TRUE)
txi.sum<-summarizeToGene(txi.kallisto,t2g,ignoreTxVersion=TRUE,countsFromAbundance="lengthScaledTPM")


write.csv(txi$abundance,asname,quote = F,row.names = T,col.names = T)
write.csv(txi.sum$abundance,asnamesum,quote = F,row.names = T,col.names = T)
write.csv(txi$counts,asname_count,quote = F,row.names = T,col.names = T)
write.csv(txi.sum$counts,asnamesum_count,quote = F,row.names = T,col.names = T)





#######利用counts数做一个deseq2比较



alltpmc <- read.csv(cname)
rownames(alltpmc) <- alltpmc$X
alltpmc$X <- NULL
CD4tpmc <- alltpmc[c(grep("CD4",colnames(alltpmc)))]
CD8tpmc <- alltpmc[c(grep("CD8",colnames(alltpmc)))]
nktpmc <- alltpmc[c(grep("_N_",colnames(alltpmc)))]
motpmc <- alltpmc[c(grep("_M_",colnames(alltpmc)))]


alltpmas <- read.csv(asname)
rownames(alltpmas) <- alltpmas$X
alltpmas$X <- NULL
CD4tpmas <- alltpmas[c(grep("CD4",colnames(alltpmas)))]
CD8tpmas <- alltpmas[c(grep("CD8",colnames(alltpmas)))]
nktpmas <- alltpmas[c(grep("_N_",colnames(alltpmas)))]
motpmas <- alltpmas[c(grep("_M_",colnames(alltpmas)))]

#错误示例这里的for条件判断语句不太对
#infoCD4 <- as.data.frame(matrix(nrow=12,ncol=0))
#infoCD4$sample <- colnames(CD4tpm)
#for (n in grep("CD4",infoCD4$sample))
#if(grep("CD4",infoCD4$sample)[n]==T)
#{
#  infoCD4$condition[n] <- "CD4"
#}

######以矩阵框的列名为依据导出条件判断矩阵
library(stringr)
infotestc <- as.data.frame(matrix(nrow=47,ncol=0))
infotestas <- as.data.frame(matrix(nrow=30,ncol=0))
infotestc$sample <- colnames(alltpmc)
infotestas$sample <- colnames(alltpmas)
infotest <- rbind(infotestc,infotestas)

#forloops <- str_detect(infotest$sample,pattern = "CD4")#创建的简单表示变量
#标注注释信息，有关细胞种类
for ( i in 1:nrow(infotest))
  
  if(str_detect(infotest$sample,pattern = "CD4")[i]==T)
  {
    infotest$condition[i] <- "CD4"
    
  }else if(str_detect(infotest$sample,pattern = "CD8")[i]==T)
  {
    infotest$condition[i] <- "CD8"
    
  }else if(str_detect(infotest$sample,pattern = "_M_")[i]==T)
  {
    infotest$condition[i] <- "MO"
    
  }else if(str_detect(infotest$sample,pattern = "_Mo_")[i]==T)
  {
    infotest$condition[i] <- "MO"
    
  }else if(str_detect(infotest$sample,pattern = "_N_")[i]==T)
  {
    infotest$condition[i] <- "NK"
    
  }else if(str_detect(infotest$sample,pattern = "_NK_")[i]==T)
  {
    infotest$condition[i] <- "NK"
  }

#标注注释信息有关组别
for ( i in 1:nrow(infotest))
  
  if(str_detect(infotest$sample,pattern = "C\\d")[i]==T)
  {
    infotest$group[i] <- "CON"
    
  }else if(str_detect(infotest$sample,pattern = "MA\\d")[i]==T)
  {
    infotest$group[i] <- "AS"
  }

write.csv(infotest,infotestname)

#对ncRNA直接用ensembl_id进行差异基因分析

cong <- read.csv(cname_count)
asg <- read.csv(asname_count)
allg <-cbind(cong,asg)
infog <- read.csv(infotestname)

allcountc <- read.csv(cname_count)

rownames(allcountc) <- allcountc$X
allcountc$X <- NULL
allcountc <- round(allcountc)
CD4countc <- allcountc[c(grep("CD4",colnames(allcountc)))]
CD8countc <- allcountc[c(grep("CD8",colnames(allcountc)))]
nkcountc <- allcountc[c(grep('_NK*_',colnames(allcountc)))]
mocountc <- allcountc[c(grep('_Mo*_',colnames(allcountc)))]
allcountas <- read.csv(asname_count)

rownames(allcountas) <- allcountas$X
allcountas$X <- NULL
allcountas <- round(allcountas)
CD4countas <- allcountas[c(grep("CD4",colnames(allcountas)))]
CD8countas <- allcountas[c(grep("CD8",colnames(allcountas)))]
nkcountas <- allcountas[c(grep("_NK*_",colnames(allcountas)))]
mocountas <- allcountas[c(grep("_Mo*_",colnames(allcountas)))]

CD4count <-cbind(CD4countc,CD4countas)
CD8count <-cbind(CD8countc,CD8countas)
nkcount <- cbind(nkcountc,nkcountas)
mocount <- cbind(mocountc,mocountas)

allcount <- cbind(allcountc,allcountas)



#写到这里是为了将分选的细胞全部隔开，然后这样或许就可以用普通的deseq流程一个个完成但是感觉工作量有点大
infog$group <- as.factor(infog$group)
infog$condition <- as.factor(infog$condition)
#######使用普通的seq流程分析的时候得到标记信息矩阵的子集
cd4infog <- subset(infog, infog$condition=="CD4")
cd8infog <- subset(infog, infog$condition=="CD8")
nkinfog <- subset(infog, infog$condition=="NK")
moinfog <- subset(infog, infog$condition=="MO")

#############使用普通的Deseq流程



all(cd4infog$sample %in% colnames(CD4count))
all(cd4infog$sample == colnames(CD4count))


dds <-  DESeqDataSetFromMatrix(countData = CD4count,colData = cd4infog,design = ~group)
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)


write.csv(diff,"deseqresult/deseq2_CD4.csv")

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diffsig,"deseqresult/deseq2_CD4_sig.csv")

up_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange < -1)

write.csv(up_DEG,"deseqresult/deseq2_CD4_sig_up.csv")
write.csv(down_DEG,"deseqresult/deseq2_CD4_sig_down.csv")


#############CD8
all(cd8infog$sample %in% colnames(CD8count))
all(cd8infog$sample == colnames(CD8count))


dds <-  DESeqDataSetFromMatrix(countData = CD8count,colData = cd8infog,design = ~group)
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)


write.csv(diff,"deseqresult/deseq2_CD8.csv")

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diffsig,"deseqresult/deseq2_CD8_sig.csv")

up_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange < -1)

write.csv(up_DEG,"deseqresult/deseq2_CD8_sig_up.csv")
write.csv(down_DEG,"deseqresult/deseq2_CD8_sig_down.csv")


###########MO




all(moinfog$sample %in% colnames(mocount))
all(moinfog$sample == colnames(mocount))


dds <-  DESeqDataSetFromMatrix(countData = mocount,colData = moinfog,design = ~group)
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)


write.csv(diff,"deseqresult/deseq2_mo.csv")

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diffsig,"deseqresult/deseq2_mo_sig.csv")

up_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange < -1)

write.csv(up_DEG,"deseqresult/deseq2_mo_sig_up.csv")
write.csv(down_DEG,"deseqresult/deseq2_mo_sig_down.csv")


###########nk



all(nkinfog$sample %in% colnames(nkcount))
all(nkinfog$sample == colnames(nkcount))


dds <-  DESeqDataSetFromMatrix(countData = nkcount,colData = nkinfog,design = ~group)
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)


write.csv(diff,"deseqresult/deseq2_nk.csv")

#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diffsig,"deseqresult/deseq2_nk_sig.csv")

up_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange < -1)

write.csv(up_DEG,"deseqresult/deseq2_nk_sig_up.csv")
write.csv(down_DEG,"deseqresult/deseq2_nk_sig_down.csv")













#cd8infog <- subset(infog,condition=="CD8")



# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.

dds$group <- factor(paste0(dds$genotype, dds$condition))
#组别信息
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

# the condition effect for genotypeIII
results(dds, contrast=c("group", "IIIB", "IIIA"))
#设定那两组比较



dds$type <- factor(paste0(dds$group,dds$condition))

design(dds) <- ~group
dds <- DESeq(dds)
resultsNames((dds))



library(DESeq2)
library(dplyr)

count <- allcount

fazdata12 <- infog
as.data.frame(fazdata12)




rownames(fazdata12) <- fazdata12$sample
all(rownames(fazdata12) %in% colnames(count))
all(rownames(fazdata12) == colnames(count))
fazdata12$condition <- as.character(fazdata12$condition)
fazdata12$group <- as.factor(fazdata12$group)
#fazdata12$patient_id <- as.character(fazdata12$patient_id)

count<- round(count)
#对于count取整如何？不知道是否正确的方法

dds <-  DESeqDataSetFromMatrix(countData = count,colData = fazdata12,design = ~group)

dds$type <- factor(paste0(dds$group,dds$condition))

design(dds) <- ~type
dds <- DESeq(dds)
resultsNames((dds))
results(dds, contrast=c("type", "AS", "CON"))
dim(dds)

#过滤
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds) 

## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)
write.csv(diff,kname2)
#Padj是P值矫正之后的数值，一般选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)

write.csv(diffsig,kname23)
up_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange > 1)
down_DEG <- subset(diffsig, padj < 0.05 & log2FoldChange < -1)

write.csv(up_DEG,kname23up)
write.csv(down_DEG,kname23down)






