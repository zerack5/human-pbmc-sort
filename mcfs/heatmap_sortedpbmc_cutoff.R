


heat2data<-read.csv("D:/data_anlysis/bulk_control_AS/bulk26/AS.casecontrol.annotated.csv", header = TRUE)



######下面是出分选细胞矩阵热图的代码，建议包成R包，可以设定几个参数变量
setwd("C:/Users/a/Desktop/表达矩阵")
"all_diffMo.csv"
"all_diffCD4.csv"
"all_diffNK.csv"
"all_diffCD8.csv"

heat2data<-read.csv("all_diffCD4.csv", header = TRUE)
heat2data<-read.csv("all_diffCD8.csv", header = TRUE)
heat2data<-read.csv("all_diffMo.csv", header = TRUE)
heat2data<-read.csv("all_diffNK.csv", header = TRUE)

colnames(heat2data)
#row.names(heat2data)<- heat2data$ext_gene
#heatdata <- heatdata[,c("ext_gene","log2FoldChange","padj","pvalue" )]
#heatdata <- heatdata[, !(names(heatdata) %in% c("ext_gene", "log2FoldChange", "padj", "pvalue"))]
#heat2data <- heat2data[,-c(1:8)]


###########@##取交集画热图
# 
# heatCD4data<-read.csv("all_diffCD4.csv", header = TRUE)
# heatCD4data <- arrange(heatCD4data,by = padj)
# heatCD4 <- heatCD4data$ext_gene[1:5000]
# 
# heatCD8data<-read.csv("all_diffCD8.csv", header = TRUE)
# heatCD8data <- arrange(heatCD8data,by = padj)
# heatCD8 <- heatCD8data$ext_gene[1:5000]
# 
# heatmodata<-read.csv("all_diffMo.csv", header = TRUE)
# heatmodata <- arrange(heatmodata,by = padj)
# heatmo <- heatmodata$ext_gene[1:5000]
# 
# heatnkdata<-read.csv("all_diffNK.csv", header = TRUE)
# heatnkdata <- arrange(heatnkdata,by = padj)
# heatnk <- heatnkdata$ext_gene[1:5000]
library(dplyr)
cutoffhold <- 1.5
threshold <- 0.01
heatCD4data<-read.csv("all_diffCD4.csv", header = TRUE)
heatCD4data <- arrange(heatCD4data,by = padj)
#heatCD4 <- heatCD4data$ext_gene[1:150]
# heatCD4data <- heatCD4data[heatCD4data$log2FoldChange > 2,]
heatCD4 <- heatCD4data$ext_gene[heatCD4data$padj < threshold]

heatCD8data<-read.csv("all_diffCD8.csv", header = TRUE)
heatCD8data <- arrange(heatCD8data,by = padj)
# heatCD8data <- heatCD8data[heatCD8data$log2FoldChange > 2,]
heatCD8 <- heatCD8data$ext_gene[heatCD8data$padj < threshold]

heatmodata<-read.csv("all_diffMo.csv", header = TRUE)
heatmodata <- arrange(heatmodata,by = padj)
# heatmodata <- heatCD4data[heatmodata$log2FoldChange > 2,]
heatmo <- heatmodata$ext_gene[heatmodata$padj < threshold]

heatnkdata<-read.csv("all_diffNK.csv", header = TRUE)
heatnkdata <- arrange(heatnkdata,by = padj)
# heatnkdata <- heatnkdata[heatnkdata$log2FoldChange > 2,]
heatnk <- heatnkdata$ext_gene[heatnkdata$padj < threshold]

cd4_cd8 <- intersect(heatCD4,heatCD8)
mo_nk <- intersect(heatmo,heatnk)
selected_gene <- intersect(cd4_cd8,mo_nk)

cutoff_CD4 <- heatCD4data$ext_gene[abs(heatCD4data$log2FoldChange) > cutoffhold]
cutoff_CD8 <- heatCD8data$ext_gene[abs(heatCD8data$log2FoldChange) > cutoffhold]
cutoff_mo <- heatmodata$ext_gene[abs(heatmodata$log2FoldChange) > cutoffhold]
cutoff_nk <- heatCD4data$ext_gene[abs(heatCD4data$log2FoldChange) > cutoffhold]

cutoff_CD4_CD8 <- intersect(cutoff_CD4,cutoff_CD8)
cutoff_mo_nk <- intersect(cutoff_mo,cutoff_nk)
cutoff_all <- intersect(cutoff_CD4_CD8,cutoff_mo_nk)

# heatCD4data <- heatCD4data[heatCD4data$ext_gene %in% selected_gene, ]
# heatCD8data <- heatCD8data[heatCD8data$ext_gene %in% selected_gene, ]
# heatmodata <- heatmodata[heatmodata$ext_gene %in% selected_gene, ]
# heatnkdata <- heatnkdata[heatnkdata$ext_gene %in% selected_gene, ]

heatCD4data <- heatCD4data[, !(names(heatCD4data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatCD8data <- heatCD8data[, !(names(heatCD8data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatmodata <- heatmodata[, !(names(heatmodata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatnkdata <- heatnkdata[, !(names(heatnkdata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]

heatCD4_CD8data <- merge(heatCD4data,heatCD8data, by = "ext_gene")
heatmo_nkdata <- merge(heatmodata,heatnkdata, by = "ext_gene")
heatalldata <- merge(heatCD4_CD8data,heatmo_nkdata, by = "ext_gene")

heatalldata <- heatalldata[heatalldata$ext_gene %in% selected_gene, ]
########为了使图像对比度更高，去掉log2foldchange小于2的值
heatalldata <- heatalldata[heatalldata$ext_gene %in% cutoff_all, ]



# 
# heatCD4data <- heatCD4data[heatCD4data$ext_gene %in% selected_gene, ]
# heatCD8data <- heatCD8data[heatCD8data$ext_gene %in% selected_gene, ]
# heatmodata <- heatmodata[heatmodata$ext_gene %in% selected_gene, ]
# heatnkdata <- heatnkdata[heatnkdata$ext_gene %in% selected_gene, ]
# 
# heatCD4data <- heatCD4data[, !(names(heatCD4data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
# heatCD8data <- heatCD8data[, !(names(heatCD8data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
# heatmodata <- heatmodata[, !(names(heatmodata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
# heatnkdata <- heatnkdata[, !(names(heatnkdata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
# 
# heatCD4_CD8data <- merge(heatCD4data,heatCD8data, by = "ext_gene")
# heatmo_nkdata <- merge(heatmodata,heatnkdata, by = "ext_gene")
# heatalldata <- merge(heatCD4_CD8data,heatmo_nkdata, by = "ext_gene")


########卸磨杀驴，将行名替换成Ext_gene并把ext_gene删掉
row.names(heatalldata) <- heatalldata$ext_gene
#heatalldata1 <- heatalldata[,!"ext_gene"]
heatalldata1 <- heatalldata[, !colnames(heatalldata) %in% "ext_gene"]


pheat123_log <- apply(heatalldata1, 2, function(x){log2(x+1)})
# pheat321 <- t(pheat123)
library(stringr)

# 创建一个空列表用于存储重复数字的映射关系
digit_mapping <- list()
counter <- 1

# 遍历每个列名
for (i in 1:length(colnames(pheat123_log))) {
  # 获取当前列名
  current_name <- colnames(pheat123_log)[i]
  
  # 提取除了"A"或"C"后面的数字以外的其他位置信息
  pattern <- "^(.*[AC])(\\d+)(.*)$"
  matches <- str_match(current_name, pattern)
  
  # 判断是否有匹配的结果
  if (nrow(matches) > 0) {
    # 获取提取的位置信息
    prefix <- matches[, 2]
    digits <- matches[, 3]
    suffix <- matches[, 4]
    
    # 判断是否有重复出现的数字
    if (!is.na(digits)) {
      # 如果数字已经在映射关系中，则直接替换为对应的整数
      if (digits %in% names(digit_mapping)) {
        colnames(pheat123_log)[i] <- paste0(prefix, digit_mapping[[digits]], suffix)
      } else {
        # 如果数字不在映射关系中，则将其替换为从1开始的连续整数，并添加到映射关系中
        colnames(pheat123_log)[i] <- paste0(prefix, counter, suffix)
        digit_mapping[[digits]] <- counter
        counter <- counter + 1
      }
    }
  }
}






##########################################################细胞分群标签
colheatdata1 <- colnames(pheat123_log)
colheatdata1<- gsub(".*_M.*", "Monocyte", colheatdata1)
colheatdata1 <- gsub(".*_N.*", "NK cell", colheatdata1)
colheatdata1 <- gsub(".*CD4.*", "CD4+Tcell", colheatdata1)
colheatdata1 <- gsub(".*CD8.*", "CD8+Tcell", colheatdata1)

annotation_col1 = data.frame(CellType = factor(colheatdata1))
ann_colors = list(
  CellType = c('Monocyte' = "#1B9E77", 'NK cell' = "#D95F02", 'CD4+Tcell' = "yellow", 'CD8+Tcell' = "pink" )
)
rownames(annotation_col1) = colnames(pheat123_log)
######################################################contro_case标签
colheatdata2 <- colnames(pheat123_log)
colheatdata2<- gsub(".*MA.*", "AS", colheatdata2)
colheatdata2 <- gsub(".*C[0-9]*.*", "Control ", colheatdata2)


annotation_col2 = data.frame(Type = factor(colheatdata2))
ann_colors = list(
  Type = c('AS' = "#D95F02", 'Control ' = "#1B9E77" )
)
rownames(annotation_col2) = colnames(pheat123_log)

#######改变合并的次序可以改变注释涂层的顺序
annotation_matrix <- cbind(annotation_col2, annotation_col1)






color_palette <- colorRampPalette(c("red", "black", "green"))(n = 500)







pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = FALSE)


pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "centroid",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette)


pheatmap::pheatmap(heatalldata1,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "centroid",
                   scale = "row",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette)



pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   scale = "row",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette)

pdf("D:/data_anlysis/bulk_control_AS/sorted_heatmap/figure/ward.D2_heatmap_test.pdf", width = 12, height = 10)
pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   scale = "row",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette,
                   annotation_colors = ann_colors,
                   annotation_col = annotation_matrix,
                   
                   )
dev.off()

#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)




######################################################################

heatdata1 <- arrange(heat2data,by = padj)
heatdata1 <- heatdata1[1:500,]
#rownames(heatdata1) <- heatdata1$ext_gene

#######根据ext_gene这一列取交集或者并集进行合并或者
merged_df <- data.frame()

base_column <- dataframes_list[[1]][[ext_gene]]

for (i in seq_along(dataframes_list)) {
  current_column <- dataframes_list[[i]][[ext_gene]]
  base_column <- intersect(base_column, current_column)
}


for (i in seq_along(dataframes_list)) {
  current_df <- dataframes_list[[i]]
  merged_df <- merge(merged_df, current_df[current_df[[common_column]] %in% base_column, ], all = TRUE)
}


heatdatat <- rbind(heatdatacd4,heatdatacd8)
heatdatamk <- rbind(heatdatamo,heatdatank)



heatdata <- heatdata1
#建议下面的直接包装成一个函数
# 创建一个逻辑向量用于标记要删除的行
to_delete <- logical(nrow(heatdata))

# 获取重复的值及其索引
duplicate_values <- heatdata[duplicated(heatdata$ext_gene), "ext_gene"]
duplicate_indexes <- which(duplicated(heatdata$ext_gene))

# 遍历每种重复值
for (i in seq_along(duplicate_values)) {
  value <- duplicate_values[i]
  indexes <- duplicate_indexes[duplicate_values == value]
  
  # 找到每种重复值中"padj"较大的行
  max_padj <- max(heatdata$padj[indexes])
  max_indexes <- indexes[heatdata$padj[indexes] == max_padj]
  
  # 将要删除的行标记为TRUE
  to_delete[max_indexes] <- TRUE
}

# 删除要删除的行
heatdata <- heatdata[!to_delete, ]

heatdata3 <- heatdata
rownames(heatdata3) <- heatdata3$ext_gene
heatdata3 <- heatdata3[,-c(1:8)]

heatdata3_log <- apply(heatdata3, 2, function(x){log2(x+1)})
pheatmap::pheatmap(heatdata3_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = FALSE)


pheatmap::pheatmap(heatdata3_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "centroid",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T)



#the agglomeration method to be used. This should be (an unambiguous abbreviation of) 
#one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).






normalize <- function(x) {
  if((max(x) - min(x)) == 0){
    return(mean(x))
  }else{
    return((x - min(x)) / (max(x) - min(x)))
  }
}
heat_map_res <- apply(results[,1:22], 2,normalize)
heat_map_res <- t(as.data.frame(heat_map_res))
mycol<-colorRampPalette(c("navy", "white", "firebrick3"))(100)
library(pheatmap)
pheatmap(heat_map_res,color = mycol,cluster_rows = F,cluster_cols = F,scale = 'none')



pheat12 <- cbind(pheat1,pheat2)
pheat123 <- cbind(pheat12,pheat3)
pheat123_log <- apply(pheat123, 2, function(x){log2(x+1)})
pheat321 <- t(pheat123)

pheatmap::pheatmap(pheat321,
                   scale = "column"
                   # color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   # kmeans_k = NA,
                   #legend = TRUE,
                   
)


pheatmap::pheatmap(pheat1236)



########################取并集画热图

library(dplyr)
# 根据padj值设置筛选条件，提取符合条件的ext_gene值


#################
cutoffhold <- 2
threshold <- 0.01
heatCD4data<-read.csv("all_diffCD4.csv", header = TRUE)
heatCD4data <- arrange(heatCD4data,by = padj)
#heatCD4 <- heatCD4data$ext_gene[1:150]
# heatCD4data <- heatCD4data[heatCD4data$log2FoldChange > 2,]
heatCD4 <- heatCD4data$ext_gene[heatCD4data$padj < threshold]

heatCD8data<-read.csv("all_diffCD8.csv", header = TRUE)
heatCD8data <- arrange(heatCD8data,by = padj)
# heatCD8data <- heatCD8data[heatCD8data$log2FoldChange > 2,]
heatCD8 <- heatCD8data$ext_gene[heatCD8data$padj < threshold]

heatmodata<-read.csv("all_diffMo.csv", header = TRUE)
heatmodata <- arrange(heatmodata,by = padj)
# heatmodata <- heatCD4data[heatmodata$log2FoldChange > 2,]
heatmo <- heatmodata$ext_gene[heatmodata$padj < threshold]

heatnkdata<-read.csv("all_diffNK.csv", header = TRUE)
heatnkdata <- arrange(heatnkdata,by = padj)
# heatnkdata <- heatnkdata[heatnkdata$log2FoldChange > 2,]
heatnk <- heatnkdata$ext_gene[heatnkdata$padj < threshold]

cd4_cd8 <- union(heatCD4,heatCD8)
mo_nk <- union(heatmo,heatnk)
selected_gene <- union(cd4_cd8,mo_nk)


cutoff_CD4 <- heatCD4data$ext_gene[abs(heatCD4data$log2FoldChange) > cutoffhold]
cutoff_CD8 <- heatCD8data$ext_gene[abs(heatCD8data$log2FoldChange) > cutoffhold]
cutoff_mo <- heatmodata$ext_gene[abs(heatmodata$log2FoldChange) > cutoffhold]
cutoff_nk <- heatCD4data$ext_gene[abs(heatCD4data$log2FoldChange) > cutoffhold]

cutoff_CD4_CD8 <- intersect(cutoff_CD4,cutoff_CD8)
cutoff_mo_nk <- intersect(cutoff_mo,cutoff_nk)
cutoff_all <- intersect(cutoff_CD4_CD8,cutoff_mo_nk)

# heatCD4data <- heatCD4data[heatCD4data$ext_gene %in% selected_gene, ]
# heatCD8data <- heatCD8data[heatCD8data$ext_gene %in% selected_gene, ]
# heatmodata <- heatmodata[heatmodata$ext_gene %in% selected_gene, ]
# heatnkdata <- heatnkdata[heatnkdata$ext_gene %in% selected_gene, ]

heatCD4data <- heatCD4data[, !(names(heatCD4data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatCD8data <- heatCD8data[, !(names(heatCD8data) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatmodata <- heatmodata[, !(names(heatmodata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]
heatnkdata <- heatnkdata[, !(names(heatnkdata) %in% c("X", "log2FoldChange", "padj", "pvalue","baseMean","lfcSE","stat"))]

heatCD4_CD8data <- merge(heatCD4data,heatCD8data, by = "ext_gene")
heatmo_nkdata <- merge(heatmodata,heatnkdata, by = "ext_gene")
heatalldata <- merge(heatCD4_CD8data,heatmo_nkdata, by = "ext_gene")

heatalldata <- heatalldata[heatalldata$ext_gene %in% selected_gene, ]
########为了使图像对比度更高，去掉log2foldchange小于2的值
heatalldata <- heatalldata[heatalldata$ext_gene %in% cutoff_all, ]


########卸磨杀驴，将行名替换成Ext_gene并把ext_gene删掉
row.names(heatalldata) <- heatalldata$ext_gene
#heatalldata1 <- heatalldata[,!"ext_gene"]
heatalldata1 <- heatalldata[, !colnames(heatalldata) %in% "ext_gene"]

# 创建自定义的红绿配色调色板
#color_palette <- colorRampPalette(c("#FF6600", "#FFFFFF", "#007FFF"))(n = 500)
color_palette <- colorRampPalette(c("red", "black", "green"))(n = 500)

pheat123_log <- apply(heatalldata1, 2, function(x){log2(x+1)})

########################################将命名修改成从一开始的连续数字
library(stringr)

# 创建一个空列表用于存储重复数字的映射关系
digit_mapping <- list()
counter <- 1

# 遍历每个列名
for (i in 1:length(colnames(pheat123_log))) {
  # 获取当前列名
  current_name <- colnames(pheat123_log)[i]
  
  # 提取除了"A"或"C"后面的数字以外的其他位置信息
  pattern <- "^(.*[AC])(\\d+)(.*)$"
  matches <- str_match(current_name, pattern)
  
  # 判断是否有匹配的结果
  if (nrow(matches) > 0) {
    # 获取提取的位置信息
    prefix <- matches[, 2]
    digits <- matches[, 3]
    suffix <- matches[, 4]
    
    # 判断是否有重复出现的数字
    if (!is.na(digits)) {
      # 如果数字已经在映射关系中，则直接替换为对应的整数
      if (digits %in% names(digit_mapping)) {
        colnames(pheat123_log)[i] <- paste0(prefix, digit_mapping[[digits]], suffix)
      } else {
        # 如果数字不在映射关系中，则将其替换为从1开始的连续整数，并添加到映射关系中
        colnames(pheat123_log)[i] <- paste0(prefix, counter, suffix)
        digit_mapping[[digits]] <- counter
        counter <- counter + 1
      }
    }
  }
}






##########################################################细胞分群标签
colheatdata1 <- colnames(pheat123_log)
colheatdata1<- gsub(".*_M.*", "Monocyte", colheatdata1)
colheatdata1 <- gsub(".*_N.*", "NK cell", colheatdata1)
colheatdata1 <- gsub(".*CD4.*", "CD4+Tcell", colheatdata1)
colheatdata1 <- gsub(".*CD8.*", "CD8+Tcell", colheatdata1)

annotation_col1 = data.frame(CellType = factor(colheatdata1))
ann_colors = list(
  CellType = c('Monocyte' = "#1B9E77", 'NK cell' = "#D95F02", 'CD4+Tcell' = "yellow", 'CD8+Tcell' = "pink" )
)
rownames(annotation_col1) = colnames(pheat123_log)
######################################################contro_case标签
colheatdata2 <- colnames(pheat123_log)
colheatdata2<- gsub(".*MA.*", "AS", colheatdata2)
colheatdata2 <- gsub(".*C[0-9]*.*", "Control ", colheatdata2)


annotation_col2 = data.frame(Type = factor(colheatdata2))
ann_colors = list(
  Type = c('AS' = "#D95F02", 'Control ' = "#1B9E77" )
)
rownames(annotation_col2) = colnames(pheat123_log)

#######改变合并的次序可以改变注释涂层的顺序
annotation_matrix <- cbind(annotation_col2, annotation_col1)






color_palette <- colorRampPalette(c("red", "black", "green"))(n = 500)





pdf("D:/data_anlysis/bulk_control_AS/sorted_heatmap/figure/complete_heatmap_test.pdf", width = 6, height = 6)

pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = FALSE)
dev.off()

pdf("D:/data_anlysis/bulk_control_AS/sorted_heatmap/figure/centroid_heatmap_test.pdf", width = 6, height = 6)
pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "centroid",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette)
dev.off()

pdf("D:/data_anlysis/bulk_control_AS/sorted_heatmap/figure/centroid_heatmap_test.pdf", width = 6, height = 6)
pheatmap::pheatmap(heatalldata1,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "centroid",
                   scale = "row",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette)

dev.off()

pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "centroid",
                   scale = "row",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette)
dev.off()

pdf("D:/data_anlysis/bulk_control_AS/sorted_heatmap/figure/ward.D2_heatmap_test.pdf", width = 6, height = 6)
pheatmap::pheatmap(pheat123_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   scale = "row",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T,
                   color = color_palette)


dev.off()







































































































































heatmap(heatalldata1, cexCol = 1,scale="row",col = greenred(100))


#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

###################################################################################设置卡padj或者log2foldchange的标准



#heatmap
library(pheatmap)

#mRNA
a<-read.table("mRNA_value.txt",header = TRUE)

annotation_col = data.frame(
  HairLossType = factor(c("Male Pattern Hair Loss","Male Pattern Hair Loss","Male Pattern Hair Loss","Male Pattern Hair Loss","Female Pattern Hair Loss"
                          ,"Female Pattern Hair Loss","Male Pattern Hair Loss","Male Pattern Hair Loss","Female Pattern Hair Loss","Female Pattern Hair Loss",
                          "Male Pattern Hair Loss","Male Pattern Hair Loss","Male Pattern Hair Loss","Male Pattern Hair Loss","Male Pattern Hair Loss",
                          "Male Pattern Hair Loss","Male Pattern Hair Loss","Male Pattern Hair Loss","Female Pattern Hair Loss","Female Pattern Hair Loss"
  ))
)

annotation_col = data.frame(
  Type = factor(c("Non-balding Area","Non-balding Area","Non-balding Area","Non-balding Area","Non-balding Area",
                  "Non-balding Area","Non-balding Area","Non-balding Area","Non-balding Area","Non-balding Area",
                  "Balding Area","Balding Area","Balding Area","Balding Area","Balding Area",
                  "Balding Area","Balding Area","Balding Area","Balding Area","Balding Area"))
)

annotation_col = data.frame(
  ScalpArea = factor(c("Occipital","Occipital","Occipital","Occipital","Occipital",
                       "Occipital","Occipital","Occipital","Occipital","Occipital",
                       "Vertex","Vertex","Vertex","Vertex","Vertex",
                       "Vertex","Vertex","Vertex","Vertex","Vertex"))
)

rownames(annotation_col) = colnames(a)
#Specify colors
ann_colors = list(
  ScalpArea = c(Occipital = "#1B9E77", Vertex = "#D95F02")
)

pheatmap(a,scale = 'row',
         show_rownames = F, border_color = "grey", 
         color = colorRampPalette(c("cyan1","black","chocolate1"))(100),
         cellwidth = 20,cellheight = 1,
         fontsize = 16, main = "Z-scored TPM Value of Differentially Expressed Genes",
         annotation_col = annotation_col, cluster_cols = F, annotation_colors = ann_colors,legend_breaks = seq(-4,3,1))









library(stringr)

# 创建一个空列表用于存储重复数字的映射关系
digit_mapping <- list()
counter <- 1

# 遍历每个列名
for (i in 1:length(colnames(pheat123_log))) {
  # 获取当前列名
  current_name <- colnames(pheat123_log)[i]
  
  # 提取除了"A"或"C"后面的数字以外的其他位置信息
  pattern <- "^(.*[AC])(\\d+)(.*)$"
  matches <- str_match(current_name, pattern)
  
  # 判断是否有匹配的结果
  if (nrow(matches) > 0) {
    # 获取提取的位置信息
    prefix <- matches[, 2]
    digits <- matches[, 3]
    suffix <- matches[, 4]
    
    # 判断是否有重复出现的数字
    if (!is.na(digits)) {
      # 如果数字已经在映射关系中，则直接替换为对应的整数
      if (digits %in% names(digit_mapping)) {
        colnames(pheat123_log)[i] <- paste0(prefix, digit_mapping[[digits]], suffix)
      } else {
        # 如果数字不在映射关系中，则将其替换为从1开始的连续整数，并添加到映射关系中
        colnames(pheat123_log)[i] <- paste0(prefix, counter, suffix)
        digit_mapping[[digits]] <- counter
        counter <- counter + 1
      }
    }
  }
}




library(stringr)

# 提取符合指定模式的数字
pattern <- "(?<=A|C)\\d+"  # 查找"A"或"C"后面紧跟着的数字
digits <- str_extract_all(colnames(pheat123_log), pattern)

# 将数字列表展平并去重
all_digits <- unlist(digits)
duplicate_digits <- all_digits[duplicated(all_digits)]

unique_duplicate_digits <- unique(duplicate_digits)
# 输出重复的数字
print(unique_duplicate_digits)












