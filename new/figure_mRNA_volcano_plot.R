
###########mRNA volcano plot
########CD4
heatdata<-read.csv("D:/as分选课题/分选课题/1/CD4_res.csv", header = TRUE)
colnames(heatdata)
heatdata <- heatdata[,c("X","log2FoldChange","padj","valpadj" )]
#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]

library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节


data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 2, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
head(data)

#概览一下选出的gene数量：
summary(data$changed)

ggplot(data,aes(log2FoldChange,-log10(padj),
                #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.3,0.9),
        legend.background = element_rect(fill='transparent'))

library(ggrepel)
#和前面一样，ifelse函数先筛选出这些基因，并形成新列：
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>2,data$X,NA)
#####data$selectedgene 从这里指定人工标注的基因
head(data)
##      Gene log2FoldChange    pvalue      padj  changed selectedgene
## 1    DOK6         0.5100 1.861e-08 0.0003053       no         <NA>
## 2    TBX5        -2.1290 5.655e-08 0.0004191 negative         TBX5
## 3 SLC32A1         0.9003 7.664e-08 0.0004191       no         <NA>
## 4  IFITM1        -1.6870 3.735e-06 0.0068090 negative       IFITM1
## 5   NUP93         0.3659 3.373e-06 0.0068090       no         <NA>
## 6 EMILIN2         1.5340 2.976e-06 0.0068090 positive      EMILIN2
ggplot(data,aes(log2FoldChange,-log10(padj),
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme(panel.grid = element_blank())+ # 去掉网格线背景
  theme(axis.line = element_line(color = "black"))+ #加上黑色的x轴与y轴
  scale_color_manual(values = c('deepskyblue','grey','tomato'))+
  scale_size_manual(values = c(2,1,2))+ 
  theme(legend.title = element_blank(),
        legend.position = c(0.5,0.9),
        legend.background = element_rect(fill='transparent'))+
  #ggrepel包内的函数给选择的基因加上文本标签:
  geom_text_repel(aes(label=selectedgene), color="black",size=4,
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black")+
  coord_equal(xlim = c(-max(abs(data$log2FoldChange)), max(abs(data$log2FoldChange))),ratio = 1)+
  theme(aspect.ratio = 1)

ggsave("CD4_m_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/as分选课题/分选课题/figure")


########CD8

heatdata<-read.csv("D:/as分选课题/分选课题/1/CD8_res.csv", header = TRUE)
colnames(heatdata)
heatdata <- heatdata[,c("X","log2FoldChange","padj","valpadj" )]
#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]

library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节


data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 2, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
head(data)

#概览一下选出的gene数量：
summary(data$changed)

ggplot(data,aes(log2FoldChange,-log10(padj),
                #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.3,0.9),
        legend.background = element_rect(fill='transparent'))

library(ggrepel)
#和前面一样，ifelse函数先筛选出这些基因，并形成新列：
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>2,data$X,NA)
#####data$selectedgene 从这里指定人工标注的基因
head(data)
##      Gene log2FoldChange    pvalue      padj  changed selectedgene
## 1    DOK6         0.5100 1.861e-08 0.0003053       no         <NA>
## 2    TBX5        -2.1290 5.655e-08 0.0004191 negative         TBX5
## 3 SLC32A1         0.9003 7.664e-08 0.0004191       no         <NA>
## 4  IFITM1        -1.6870 3.735e-06 0.0068090 negative       IFITM1
## 5   NUP93         0.3659 3.373e-06 0.0068090       no         <NA>
## 6 EMILIN2         1.5340 2.976e-06 0.0068090 positive      EMILIN2
ggplot(data,aes(log2FoldChange,-log10(padj),
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  scale_color_manual(values = c('deepskyblue','grey','tomato'))+
  scale_size_manual(values = c(2,1,2))+ 
  theme(legend.title = element_blank(),
        legend.position = c(0.5,0.9),
        legend.background = element_rect(fill='transparent'))+
  #ggrepel包内的函数给选择的基因加上文本标签:
  geom_text_repel(aes(label=selectedgene), color="black",size=4,
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 

ggsave("CD8_m_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/as分选课题/分选课题/figure")

 
####################Mo
heatdata<-read.csv("D:/as分选课题/分选课题/1/res_Mo.csv", header = TRUE)
colnames(heatdata)
heatdata$padj <- p.adjust(heatdata$pvalue, method = "fdr")
heatdata <- heatdata[,c("X","log2FoldChange","padj","valpadj" )]
#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]

library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节


data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 2, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
head(data)

#概览一下选出的gene数量：
summary(data$changed)

ggplot(data,aes(log2FoldChange,-log10(padj),
                #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.3,0.9),
        legend.background = element_rect(fill='transparent'))

library(ggrepel)
#和前面一样，ifelse函数先筛选出这些基因，并形成新列：
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>2,data$X,NA)
#####data$selectedgene 从这里指定人工标注的基因
head(data)
##      Gene log2FoldChange    pvalue      padj  changed selectedgene
## 1    DOK6         0.5100 1.861e-08 0.0003053       no         <NA>
## 2    TBX5        -2.1290 5.655e-08 0.0004191 negative         TBX5
## 3 SLC32A1         0.9003 7.664e-08 0.0004191       no         <NA>
## 4  IFITM1        -1.6870 3.735e-06 0.0068090 negative       IFITM1
## 5   NUP93         0.3659 3.373e-06 0.0068090       no         <NA>
## 6 EMILIN2         1.5340 2.976e-06 0.0068090 positive      EMILIN2
ggplot(data,aes(log2FoldChange,-log10(padj),
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  scale_color_manual(values = c('deepskyblue','grey','tomato'))+
  scale_size_manual(values = c(2,1,2))+ 
  theme(legend.title = element_blank(),
        legend.position = c(0.5,0.9),
        legend.background = element_rect(fill='transparent'))+
  #ggrepel包内的函数给选择的基因加上文本标签:
  geom_text_repel(aes(label=selectedgene), color="black",size=4,
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 

ggsave("Mo_m_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/as分选课题/分选课题/figure")

###########################NK
heatdata<-read.csv("D:/as分选课题/分选课题/1/res_NK.csv", header = TRUE)
colnames(heatdata)
heatdata <- heatdata[,c("X","log2FoldChange","padj","valpadj" )]
#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]

library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节


data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.001 & abs(data$log2FoldChange) > 2, 
                              ifelse(data$log2FoldChange > 1,'positive','negative'),'no'))
head(data)

#概览一下选出的gene数量：
summary(data$changed)

ggplot(data,aes(log2FoldChange,-log10(padj),
                #分别给正负显著变化的基因在图中根据颜色、大小标注出来：
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  scale_color_manual(values = c('blue','grey','red'))+
  scale_size_manual(values = c(2,1,2))+
  theme(legend.title = element_blank(), #图例的设置参数
        legend.position = c(0.3,0.9),
        legend.background = element_rect(fill='transparent'))

library(ggrepel)
#和前面一样，ifelse函数先筛选出这些基因，并形成新列：
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>2,data$X,NA)
#####data$selectedgene 从这里指定人工标注的基因
head(data)
##      Gene log2FoldChange    pvalue      padj  changed selectedgene
## 1    DOK6         0.5100 1.861e-08 0.0003053       no         <NA>
## 2    TBX5        -2.1290 5.655e-08 0.0004191 negative         TBX5
## 3 SLC32A1         0.9003 7.664e-08 0.0004191       no         <NA>
## 4  IFITM1        -1.6870 3.735e-06 0.0068090 negative       IFITM1
## 5   NUP93         0.3659 3.373e-06 0.0068090       no         <NA>
## 6 EMILIN2         1.5340 2.976e-06 0.0068090 positive      EMILIN2
ggplot(data,aes(log2FoldChange,-log10(padj),
                color=factor(changed),
                size=factor(changed)))+  
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)"))+
  theme_grey(base_size = 15)+
  scale_color_manual(values = c('deepskyblue','grey','tomato'))+
  scale_size_manual(values = c(2,1,2))+ 
  theme(legend.title = element_blank(),
        legend.position = c(0.5,0.9),
        legend.background = element_rect(fill='transparent'))+
  #ggrepel包内的函数给选择的基因加上文本标签:
  geom_text_repel(aes(label=selectedgene), color="black",size=4,
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black") 

ggsave("nk_m_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/as分选课题/分选课题/figure")


