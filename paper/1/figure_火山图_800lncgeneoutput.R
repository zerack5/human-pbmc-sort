

path0 <- "D:/data_anlysis/asncRNA"
setwd(path0)


file1 <- "1phase_il17acount.csv"
file2 <- "2phase_il17acount.csv"
file3 <- "3phase_il17acount.csv"

de_cd4 <- "deseqresult/deseq2_CD4.csv"
de_cd4_sig <- "deseqresult/deseq2_CD4_sig.csv"
de_cd4_sig_up <- "deseqresult/deseq2_CD4_sig_up.csv"
de_cd4_sig_down <- "deseqresult/deseq2_CD4_sig_down.csv"

de_cd8 <- "deseqresult/deseq2_CD8.csv"
de_cd8_sig <- "deseqresult/deseq2_CD8_sig.csv"
de_cd8_sig_up <- "deseqresult/deseq2_CD8_sig_up.csv"
de_cd8_sig_down <- "deseqresult/deseq2_CD8_sig_down.csv"

de_mo <- "deseqresult/deseq2_mo.csv"
de_mo_sig <- "deseqresult/deseq2_mo_sig.csv"
de_mo_sig_up <- "deseqresult/deseq2_mo_sig_up.csv"
de_mo_sig_down <- "deseqresult/deseq2_mo_sig_down.csv"

de_nk <- "deseqresult/deseq2_nk.csv"
de_nk_sig <- "deseqresult/deseq2_nk_sig.csv"
de_nk_sig_up <- "deseqresult/deseq2_nk_sig_up.csv"
de_nk_sig_down <- "deseqresult/deseq2_nk_sig_down.csv"




diff<-read.csv(file=de_cd4,sep = ",")

library(dplyr)
diff <- arrange(diff,padj)

colnames(diff)[1] <- "ensembl_gene_id"

diff$log2FoldChange <- (-1)*(diff$log2FoldChange)

#########volcano plot of deseq result of CD4

heat1data<-diff
colnames(heat1data)
heatdata <- heat1data[,c("external_gene_name","log2FoldChange","pvalue","padj" )]
library(EnhancedVolcano)
EnhancedVolcano(heat1data,
                lab = heat1data$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'difference between phase1 and phase2',
                pCutoff = 10e-4,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)

#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]
library(ggrepel)
library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节

data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 1, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
data$changed <- as.character(data$changed)
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>1,data$external_gene_name,NA)
head(data)
#概览一下选出的gene数量：
summary(data$changed)
ggplot(data=data,aes(x =log2FoldChange,y = -log10(padj),
                     #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                     color=factor(changed),
                     size=factor(changed)))+  
  geom_point(alpha=0.8 , size =1 ,shape =19)+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  ggtitle("CD4 leftcontrol rightas")+
  
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(panel.background = element_rect(fill = "gray90", color = "white"),#改变图表面板（panel）的背景颜色和样式
        panel.grid.major = element_blank(),#s删除散点图的主要网格线
        panel.grid.minor = element_blank(),#删除散点图中的次要网格线
        legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.1,0.9),
        legend.background = element_rect(fill='transparent'))+
  geom_hline(yintercept=-log10(0.001),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  geom_text_repel(aes(label=selectedgene), color="black",size=3,#为显著的基因基于算法不重合地添加genename
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 

ggsave("CD4_lnc_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/data_anlysis/asncRNA/figure_test/lnc/")
#########volcano plot of deseq result of CD8

diff<-read.csv(file=de_cd8,sep = ",")

library(dplyr)
diff <- arrange(diff,padj)
diff$log2FoldChange <- (-1)*(diff$log2FoldChange)

colnames(diff)[1] <- "ensembl_gene_id"
heat1data<-diff
colnames(heat1data)
heatdata <- heat1data[,c("external_gene_name","log2FoldChange","pvalue","padj" )]
library(EnhancedVolcano)
EnhancedVolcano(heat1data,
                lab = heat1data$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'difference between case and control',
                pCutoff = 10e-4,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)

#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]
library(ggrepel)
library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节

data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 1, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
data$changed <- as.character(data$changed)
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>1,data$external_gene_name,NA)
head(data)
#概览一下选出的gene数量：
summary(data$changed)
ggplot(data=data,aes(x =log2FoldChange,y = -log10(padj),
                     #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                     color=factor(changed),
                     size=factor(changed)))+  
  geom_point(alpha=0.8 , size =1 ,shape =19)+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  ggtitle("CD8 leftcontrol rightas")+
  
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(panel.background = element_rect(fill = "gray90", color = "white"),#改变图表面板（panel）的背景颜色和样式
        panel.grid.major = element_blank(),#s删除散点图的主要网格线
        panel.grid.minor = element_blank(),#删除散点图中的次要网格线
        legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.1,0.9),
        legend.background = element_rect(fill='transparent'))+
  geom_hline(yintercept=-log10(0.001),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  geom_text_repel(aes(label=selectedgene), color="black",size=3,#为显著的基因基于算法不重合地添加genename
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 

ggsave("CD8_lnc_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/data_anlysis/asncRNA/figure_test/lnc/")


#########volcano plot of deseq result of mo

diff<-read.csv(file=de_mo,sep = ",")

library(dplyr)
diff <- arrange(diff,padj)
diff$log2FoldChange <- (-1)*(diff$log2FoldChange)

colnames(diff)[1] <- "ensembl_gene_id"
heat1data<-diff
colnames(heat1data)
heatdata <- heat1data[,c("external_gene_name","log2FoldChange","pvalue","padj" )]
library(EnhancedVolcano)
EnhancedVolcano(heat1data,
                lab = heat1data$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'difference between case and control',
                pCutoff = 10e-4,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)

#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]
library(ggrepel)
library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节

data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 1, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
data$changed <- as.character(data$changed)
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>1,data$external_gene_name,NA)
head(data)
#概览一下选出的gene数量：
summary(data$changed)
ggplot(data=data,aes(x =log2FoldChange,y = -log10(padj),
                     #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                     color=factor(changed),
                     size=factor(changed)))+  
  geom_point(alpha=0.8 , size =1 ,shape =19)+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  ggtitle("MO leftcontrol rightas")+
  
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(panel.background = element_rect(fill = "gray90", color = "white"),#改变图表面板（panel）的背景颜色和样式
        panel.grid.major = element_blank(),#s删除散点图的主要网格线
        panel.grid.minor = element_blank(),#删除散点图中的次要网格线
        legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.1,0.9),
        legend.background = element_rect(fill='transparent'))+
  geom_hline(yintercept=-log10(0.001),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  geom_text_repel(aes(label=selectedgene), color="black",size=3,#为显著的基因基于算法不重合地添加genename
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 
ggsave("mo_lnc_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/data_anlysis/asncRNA/figure_test/lnc/")

#########volcano plot of deseq result of nk

diff<-read.csv(file=de_nk,sep = ",")

library(dplyr)
diff <- arrange(diff,padj)
diff$log2FoldChange <- (-1)*(diff$log2FoldChange)

colnames(diff)[1] <- "ensembl_gene_id"
heat1data<-diff
colnames(heat1data)
heatdata <- heat1data[,c("external_gene_name","log2FoldChange","pvalue","padj" )]
library(EnhancedVolcano)
EnhancedVolcano(heat1data,
                lab = heat1data$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'difference between case and control',
                pCutoff = 10e-4,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)

#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]
library(ggrepel)
library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节

data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 1, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
data$changed <- as.character(data$changed)
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>1,data$external_gene_name,NA)
head(data)
#概览一下选出的gene数量：
summary(data$changed)
ggplot(data=data,aes(x =log2FoldChange,y = -log10(padj),
                     #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                     color=factor(changed),
                     size=factor(changed)))+  
  geom_point(alpha=0.8 , size =1 ,shape =19)+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  ggtitle("Nk leftcontrol rightas")+
  
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(panel.background = element_rect(fill = "gray90", color = "white"),#改变图表面板（panel）的背景颜色和样式
        panel.grid.major = element_blank(),#s删除散点图的主要网格线
        panel.grid.minor = element_blank(),#删除散点图中的次要网格线
        legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.1,0.9),
        legend.background = element_rect(fill='transparent'))+
  geom_hline(yintercept=-log10(0.001),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  geom_text_repel(aes(label=selectedgene), color="black",size=3,#为显著的基因基于算法不重合地添加genename
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 


ggsave("nk_lnc_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/data_anlysis/asncRNA/figure_test/lnc/")




#library(biomaRt)
#listMarts()
#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                         dataset = "hsapiens_gene_ensembl",
#                         host = 'ensembl.org')
#t5g <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","description","chromosome_name"), mart = mart)
#saveRDS(t5g,"D:/Rdata/biomart/t5g.rdata")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

diff<-read.csv(file=de_cd4,sep = ",")

library(dplyr)
diff <- arrange(diff,padj)

colnames(diff)[1] <- "ensembl_gene_id"


t5g <- readRDS("D:/Rdata/biomart/t5g.rdata")


tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]

#下面是biomart拉出较多选项时留出的草稿
#tx2g<- t2g
#tx2g$ensembl_transcript_id <- NULL
tx2g <- unique(tx2g)


dataf <- merge(diff,tx2g,by="ensembl_gene_id")

dataf <- arrange(dataf,desc(log2FoldChange))

#tx2g$ensembl_gene_id
#diff$ensembl_gene_id




#划定用来做富集的基因数量
#datak <- subset(dataf,dataf$pvalue<0.05)
#datak <- subset(dataf,dataf$padj<0.05)
#datak <- arrange(datak,padj)



#########采用NORCE包



library(NoRCE)

cd4gene_ncRNA <- dataf$ensembl_gene_id
cd4gene_ncRNA <- cd4gene_ncRNA[0:800]

setParameters("searchRegion", "exon")

#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = cd4gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')

diff<-read.csv(file=de_cd4_sig_up,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
cd4gene_ncRNA <- dataf$ensembl_gene_id
cd4gene_ncRNA <- cd4gene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = cd4gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_cd4_sig_up.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = cd4gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_cd4_sig_up.rdata")




diff<-read.csv(file=de_cd4_sig_down,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
cd4gene_ncRNA <- dataf$ensembl_gene_id
cd4gene_ncRNA <- cd4gene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = cd4gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_cd4_sig_down.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = cd4gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_cd4_sig_down.rdata")


##########CD8

diff<-read.csv(file=de_cd8,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
cd8gene_ncRNA <- dataf$ensembl_gene_id
cd8gene_ncRNA <- cd8gene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = cd8gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')

#up应该是control里面高
diff<-read.csv(file=de_cd8_sig_up,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
cd8gene_ncRNA <- dataf$ensembl_gene_id
cd8gene_ncRNA <- cd8gene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = cd8gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_cd8_sig_up.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = cd8gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_cd8_sig_up.rdata")



diff<-read.csv(file=de_cd8_sig_down,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
cd8gene_ncRNA <- dataf$ensembl_gene_id
cd8gene_ncRNA <- cd8gene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = cd8gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_cd8_sig_down.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = cd8gene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_cd8_sig_down.rdata")



#########MO

diff<-read.csv(file=de_mo,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
mogene_ncRNA <- dataf$ensembl_gene_id
mogene_ncRNA <- mogene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = mogene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')

######up意味着在contro中高表达
diff<-read.csv(file=de_mo_sig_up,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
mogene_ncRNA <- dataf$ensembl_gene_id
mogene_ncRNA <- mogene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = mogene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_mo_sig_up.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = mogene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_mo_sig_up.rdata")



#####down意味着在as中高表达
diff<-read.csv(file=de_mo_sig_down,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
mogene_ncRNA <- dataf$ensembl_gene_id
mogene_ncRNA <- mogene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = mogene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_mo_sig_down.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = mogene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_mo_sig_down.rdata")






#######NK
diff<-read.csv(file=de_nk,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
nkgene_ncRNA <- dataf$ensembl_gene_id
nkgene_ncRNA <- nkgene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = nkgene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')




#####de
diff<-read.csv(file=de_nk_sig_up,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
nkgene_ncRNA <- dataf$ensembl_gene_id
nkgene_ncRNA <- nkgene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = nkgene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_nk_sig_up.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = nkgene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_nk_sig_up.rdata")




#####down
diff<-read.csv(file=de_nk_sig_down,sep = ",")
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
library(NoRCE)
nkgene_ncRNA <- dataf$ensembl_gene_id
nkgene_ncRNA <- nkgene_ncRNA[0:800]
setParameters("searchRegion", "exon")
#NoRCE automatically consider only the exon of the genes
ncGO<-geneGOEnricher(gene = nkgene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene')
saveRDS(ncGO,"ncgoresult/de_nk_sig_down.rdata")
setParameters("pathwayType", "kegg")
ncKEGG <- genePathwayEnricher(gene = nkgene_ncRNA, org_assembly='hg38', near=TRUE, genetype = 'Ensembl_gene' )
saveRDS(ncKEGG,"nckeggresult/de_nk_sig_down.rdata")

ncCD4c_up <- readRDS("ncgoresult/de_cd4_sig_up.rdata")
ncCD4c_down <- readRDS("ncgoresult/de_cd4_sig_down.rdata")
ncCD8c_up <- readRDS("ncgoresult/de_cd8_sig_up.rdata")
ncCD8c_down <- readRDS("ncgoresult/de_cd8_sig_down.rdata")
ncmoc_up <- readRDS("ncgoresult/de_mo_sig_up.rdata")
ncmoc_down <- readRDS("ncgoresult/de_mo_sig_down.rdata")
ncnkc_up <- readRDS("ncgoresult/de_nk_sig_up.rdata")
ncnkc_down <- readRDS("ncgoresult/de_nk_sig_down.rdata")


#drawDotPlot(mrnaObject = ncCD4c_up, type = "pAdjust", 25)
#getGoDag(mrnaObject = ncCD4c_up, type = "pAdjust", n = 5,filename = 'getgodag/ncCD4c_up.png',imageFormat = "png", p_range = seq(0, 0.001, by = 0.001))

#ncCD4c_up <- readRDS("ncgoresult/1000基因输入/de_cd4_sig_up.rdata")

#######kegg读这行
#ncCD4c_up <- readRDS("nckeggresult/de_cd4_sig_up.rdata")
#ncCD4c_down <- readRDS("nckeggresult/de_cd4_sig_down.rdata")
#ncCD8c_up <- readRDS("nckeggresult/de_cd8_sig_up.rdata")
#ncCD8c_down <- readRDS("nckeggresult/de_cd8_sig_down.rdata")
#ncmoc_up <- readRDS("nckeggresult/de_mo_sig_up.rdata")
#ncmoc_down <- readRDS("nckeggresult/de_mo_sig_down.rdata")
#ncnkc_up <- readRDS("nckeggresult/de_nk_sig_up.rdata")
#ncnkc_down <- readRDS("nckeggresult/de_nk_sig_down.rdata")








tncCD4c_up <- ncCD4c_up@Term
tncCD4c_down <- ncCD4c_down@Term
tncCD8c_up <- ncCD8c_up@Term
tncCD8c_down <- ncCD8c_down@Term
tncmoc_up <- ncmoc_up@Term
tncmoc_down <- ncmoc_down@Term
tncnkc_up <- ncnkc_up@Term
tncnkc_down <- ncnkc_down@Term

CD4up_CD8up <- intersect(tncCD4c_up, tncCD8c_up)
CD4down_CD8down <- intersect(tncCD4c_down,tncCD8c_down)

moup_nkup <- intersect(tncmoc_up,tncnkc_up)
modown_nkdown <- intersect(tncmoc_down,tncnkc_down)

CD4up_CD8up_moup_nkup <- intersect(CD4up_CD8up,moup_nkup)
CD4down_CD8down_modown_nkdown <- intersect(CD4down_CD8down,modown_nkdown)
write.table(CD4up_CD8up_moup_nkup,file = "ncgoresult/CD4up_CD8up_moup_nkup03.txt",sep = "\t",row.names = F,col.names = F)
write.table(CD4down_CD8down_modown_nkdown,file = "ncgoresult/CD4down_CD8down_modown_nkdown03.txt",sep = "\t",row.names = F,col.names = F)

write.table(tncCD8c_up,file = "ncgoresult/tncCD8c_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncCD8c_down,file = "ncgoresult/tncCD8c_down.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncCD4c_up,file = "ncgoresult/tncCD4c_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncCD4c_down,file = "ncgoresult/tncCD4c_down.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncmoc_up,file = "ncgoresult/tncmoc_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncmoc_down,file = "ncgoresult/tncmoc_down.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncnkc_up,file = "ncgoresult/tncnkc_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncnkc_down,file = "ncgoresult/tncnkc_down.txt",sep = "\t",row.names = F,col.names = F)





#kegg读这里
write.table(CD4up_CD8up_moup_nkup,file = "nckeggresult/CD4up_CD8up_moup_nkup03.txt",sep = "\t",row.names = F,col.names = F)
write.table(CD4down_CD8down_modown_nkdown,file = "nckeggresult/CD4down_CD8down_modown_nkdown03.txt",sep = "\t",row.names = F,col.names = F)
write.table(CD4up_CD8up,file = "nckeggresult/CD4up_CD8up03.txt",sep = "\t",row.names = F,col.names = F)
write.table(CD4down_CD8down,file = "nckeggresult/CD4down_CD8down03.txt",sep = "\t",row.names = F,col.names = F)

write.table(tncCD8c_up,file = "nckeggresult/tncCD8c_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncCD8c_down,file = "nckeggresult/tncCD8c_down.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncCD4c_up,file = "nckeggresult/tncCD4c_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncCD4c_down,file = "nckeggresult/tncCD4c_down.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncmoc_up,file = "nckeggresult/tncmoc_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncmoc_down,file = "nckeggresult/tncmoc_down.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncnkc_up,file = "nckeggresult/tncnkc_up.txt",sep = "\t",row.names = F,col.names = F)
write.table(tncnkc_down,file = "nckeggresult/tncnkc_down.txt",sep = "\t",row.names = F,col.names = F)

createNetwork(mrnaObject = ncCD8c_down, type = 'pvalue', n = 2, isNonCode = TRUE)
#####报错信息提示target edge = 0,说明没有连到任何的点，考虑是输入的genelist太短

setParameters("pathwayType","reactome")






ncGO2 <- genePathwayEnricher(gene =  )






#gene_fc <- dataf$log2FoldChange #把foldchange按照从大到小提取出来
#head(gene_fc)
#names(gene_fc) <- dataf$entrezgene_id#给上面提取的foldchange加上对应上ENTREZID
#head(gene_fc)
#names(gene_fc)[1] <- 28452
#head(gene_fc)
#gene_fc <- na.omit(gene_fc)
#gene_fc

#genefc <- gene_fc[-1]
#head(genefc)

datah <-dataf[,c("entrezgene_id","log2FoldChange")]
head(datah)
datae <- na.omit(datah)#去除NA值因为有些当时基因名用biomart包转化后没有对上entreze_id
gene_fc <- datae$log2FoldChange
names(gene_fc) <- as.character(datae$entrezgene_id)
gene_fc = sort(gene_fc, decreasing = TRUE)#排一下序输入


#BiocManager::install("pathview")

#BiocManager::install("clusterProfiler")
#BiocManager::install("GSEABase")

#install.packages("stats")
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(DOSE)
library(GSEABase)
library(stats)

kk2 <- gseKEGG(geneList     = gene_fc,
               organism     = 'hsa',
               #nPerm        = 10000,
               pvalueCutoff = 0.05)
kk2@result <- arrange(kk2@result, desc(kk2@result$enrichmentScore))

View(kk2@result) 

gseaplot2(kk2,1:3,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值


ggplot(kk2@result, aes(x = Description, y = enrichmentScore, color = enrichmentScore > 0)) + 
  geom_line(size = 1) +
  geom_vline(xintercept = which(kk2@result$leading_edge == 1), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("blue", "red"), labels = c("Non-significant", "Significant")) +
  xlab("Ranked list") +
  ylab("Enrichment score") +
  ggtitle("GSEA Enrichment Plot") +
  theme_bw()


res <- gseGO(
  gene_fc,    # 根据logFC排序的基因集
  ont = "ALL",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
  OrgDb = org.Hs.eg.db,    # 使用人的OrgDb
  keyType = "ENSEMBL",    # 基因id类型
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",    # p值校正方法
)

head(res,12)





