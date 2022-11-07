heatdata<-read.csv("C:/Users/a/Desktop/as分选课题/分选课题/1/res_Mo.csv", header = TRUE)
colnames(heatdata)
heatdata <- heatdata[,c("X","log2FoldChange","valpadj","padj" )]
#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]

library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节


data <- heatdata
data$changed <- factor(ifelse(data$padj < 0.000001 & abs(data$log2FoldChange) > 3, 
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
data$selectedgene <- ifelse(data$padj<0.0000001 & abs(data$log2FoldChange)>3,data$X,NA)
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

