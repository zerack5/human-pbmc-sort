


heatdata<-read.csv("D:/data_anlysis/bulk_control_AS/bulk26/AS.casecontrol.annotated.csv", header = TRUE)
colnames(heatdata)
heatdata <- heatdata[,c("ext_gene","log2FoldChange","padj","pvalue" )]
#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]

library(ggplot2)
ggplot(heatdata,aes(log2FoldChange,-log10(padj)))+ 
  geom_point()+
  labs(x=expression(Log[2]*" Fold Change"),
       y=expression(-Log[10]*" (padj)")) #修改坐标轴命名细节

data <- heatdata[!is.na(heatdata$log2FoldChange), ]
data$changed <- factor(ifelse(data$padj < 0.01 & abs(data$log2FoldChange) > 2, 
                              ifelse(data$log2FoldChange > 1,'Up regulated','Down regulated'),'Not significience'))
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
        legend.position = c(0.3,0.9),###position改一下
        legend.background = element_rect(fill='transparent'))
#######legend单独放右边
library(ggrepel)
#和前面一样，ifelse函数先筛选出这些基因，并形成新列：

###阈值调一下
data$selectedgene <- ifelse(data$padj<0.001 & abs(data$log2FoldChange)>2,data$ext_gene,NA)
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
       y=expression(-Log[10]*" (padj)"),
       family = "sans", face = "bold"
       )+
  theme(panel.grid = element_blank(),# 去掉网格线背景
        axis.line = element_line(color = "black"),#加上黑色的x轴与y轴
        text = element_text(family = "sans", face = "bold"))+  #Arial字体并且加粗
  scale_color_manual(values = c('deepskyblue','grey','tomato'))+
  scale_size_manual(values = c(2,1,2))+ 
  theme(legend.title = element_blank(),
        legend.justification = "right",  # 将图例放置在右边
        #legend.position = c(0.5,0.9),
        panel.background = element_rect(fill = "white"))+  # 设置白色背景
  #ggrepel包内的函数给选择的基因加上文本标签:
  geom_text_repel(aes(label=selectedgene), color="black",size=4,
                  box.padding=unit(0.5, "lines"), 
                  point.padding=NA, 
                  segment.colour = "black")+
  coord_equal(xlim = c(-max(abs(data$log2FoldChange))*0.5, max(abs(data$log2FoldChange))*0.5),ratio = 1)+
  theme(aspect.ratio = 1,
        text = element_text(family = "sans", face = "bold"))
  

ggsave("PBMC26_m_leftcontrol_rightas.tiff",device = "tiff",width =7,height =7,dpi = 300,path = "D:/data_anlysis/bulk_control_AS/bulk26/figure")


############热图





#############cibersort_analysis
heatdata<-read.csv("D:/data_anlysis/bulk_control_AS/bulk26/AS.casecontrol.annotated.csv", header = TRUE)
mat <- heatdata[,-c(2:8)]
matrix <- unique(mat)

matrix_mean=aggregate(.~ext_gene,mean,data=matrix)
write.table(matrix_mean,file = "casecontrol_matforcibersort.txt",quote = F,sep = "\t",row.names = F)

#duplicate_rows <- duplicated(matrix$ext_gene)
#for(i in 1:length(duplicate_rows))
#if(matrix$ext_gene[i] == T)
#  print("有重复")




setwd("D:/data_anlysis/bulk_control_AS/bulk26/cibersort")



library(org.Hs.eg.db)

library("preprocessCore")

source('Cibersort.R')
result1 <- CIBERSORT("cibersortsigmatrix/486hcK999sortedcell.txt","casecontrol_matforcibersort.txt",perm = 1000, QN = T)  #perm 
write.csv(result1,file = "LM22_casecontrol_cibersortresult.csv")


source('Cibersort.R')
result2 <- CIBERSORT("cibersortsigmatrix/LM22.txt","casecontrol_matforcibersort.txt",perm = 1000, QN = T)  #perm 
write.csv(result2,file = "Cell486_casecontrol_cibersortresult.csv")

#将矩阵根据列名中是否含有c一分为2
# 指定要查找的字符

df <- result1
search_char <- "example"
# 找到列名含有指定字符的列
cols_with_char <- grep(search_char, colnames(df), value = TRUE)
# 找到列名不含指定字符的列
cols_without_char <- colnames(df)[!grepl(search_char, colnames(df))]
# 根据列名创建子集
result1_with_char <- df[, cols_with_char]
result1_without_char <- df[, cols_without_char]

df <- result2
search_char <- "example"
# 找到列名含有指定字符的列
cols_with_char <- grep(search_char, colnames(df), value = TRUE)
# 找到列名不含指定字符的列
cols_without_char <- colnames(df)[!grepl(search_char, colnames(df))]
# 根据列名创建子集
result2_with_char <- df[, cols_with_char]
result2_without_char <- df[, cols_without_char]


##########486版细胞比例

res <- data.frame(result1[,1:27])%>%
  mutate(group = c(rep('control',26),rep('as',26)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:28],
               names_to = "cell.type",
               values_to = 'value')
head(res,n=42)
res$group <- factor(res$group, levels = c('control', 'as'))
############对于kruskal.test后的结果进行筛选，只保留P值符合要求的组别
kruskal_result <- res %>%
  group_by(cell.type) %>%
  summarise(p_value = kruskal.test(value ~ group)$p.value)

res$pvalue <- NA
for(i in unique(res$cell.type)){
  the_dt <- res[which(res$cell.type == i),]
  res$pvalue[which(res$cell.type == i)] <- kruskal.test(data = the_dt, value ~ group)$p.value
}

filtered_res <- filter(res, pvalue <= 0.05)  #展示的图形数量阈值设定
############################################################################

ggplot(filtered_res,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  #scale_fill_nejm()+
  #scale_fill_nejm()是ggplot2中一个用于调整填充颜色的函数，它用于设置图形中填充颜色的默认值，使其与New England Journal of Medicine（NEJM）杂志的风格一致。
  
  scale_fill_manual(values = c("control" = "green",
                               "as" = "red" )) +
  #要调整箱型图分组的颜色，可以在 scale_fill_nejm() 中使用 manual 选项来指定颜色
  
  
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")+
  #使用 stat_compare_means() 函数在图表上添加统计学检验的结果，其中 group 是分组变量，使用 "kruskal.test" 方法计算组间差异。检验结果使用显著性水平标注在图表上。
  
  ########将表格当中的用药前中后三种组别的差异分别做三次检验来体现出来
  geom_signif(comparisons = list(c("control", "as")),
              y_position = c(0.9, 0.8),
              annotations = c("**", "*"),
              tip_length = 0.02)

geom_signif(comparisons = list(c("control", "as")), map_signif_level = TRUE, y_position = 0.9) 

#############LM22版细胞比例


res <- data.frame(result2[,1:22])%>%
  mutate(group = c(rep('control',26),rep('as',26)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value')
head(res,n=42)
res$group <- factor(res$group, levels = c('control', 'as'))
############对于kruskal.test后的结果进行筛选，只保留P值符合要求的组别
kruskal_result <- res %>%
  group_by(cell.type) %>%
  summarise(p_value = kruskal.test(value ~ group)$p.value)

res$pvalue <- NA
for(i in unique(res$cell.type)){
  the_dt <- res[which(res$cell.type == i),]
  res$pvalue[which(res$cell.type == i)] <- kruskal.test(data = the_dt, value ~ group)$p.value
}

filtered_res <- filter(res, pvalue < 0.05)  ####展示的图形阈值设定
############################################################################

ggplot(filtered_res,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  #scale_fill_nejm()+
  #scale_fill_nejm()是ggplot2中一个用于调整填充颜色的函数，它用于设置图形中填充颜色的默认值，使其与New England Journal of Medicine（NEJM）杂志的风格一致。
  
  scale_fill_manual(values = c("control" = "green",
                               "as" = "red" )) +
  #要调整箱型图分组的颜色，可以在 scale_fill_nejm() 中使用 manual 选项来指定颜色
  
  
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")+
  #使用 stat_compare_means() 函数在图表上添加统计学检验的结果，其中 group 是分组变量，使用 "kruskal.test" 方法计算组间差异。检验结果使用显著性水平标注在图表上。
  
  ########将表格当中的用药前中后三种组别的差异分别做三次检验来体现出来
  geom_signif(comparisons = list(c("control", "as")),
              y_position = c(0.9, 0.8),
              annotations = c("**", "*"),
              tip_length = 0.02)

geom_signif(comparisons = list(c("control", "as")), map_signif_level = TRUE, y_position = 0.9) 


