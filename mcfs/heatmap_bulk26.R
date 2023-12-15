


heat2data<-read.csv("D:/data_anlysis/bulk_control_AS/bulk26/AS.casecontrol.annotated.csv", header = TRUE)



######下面是出分选细胞矩阵热图的代码，建议包成R包，可以设定几个参数变量
#setwd("C:/Users/a/Desktop/表达矩阵")
#"all_diffMo.csv"
#"all_diffCD4.csv"
#"all_diffNK.csv"
#"all_diffCD8.csv"

#heat2data<-read.csv("all_diffCD4.csv", header = TRUE)
#heat2data<-read.csv("all_diffCD8.csv", header = TRUE)
#heat2data<-read.csv("all_diffMo.csv", header = TRUE)
#heat2data<-read.csv("all_diffNK.csv", header = TRUE)

colnames(heat2data)
#row.names(heat2data)<- heat2data$ext_gene
#heatdata <- heatdata[,c("ext_gene","log2FoldChange","padj","pvalue" )]
#heatdata <- heatdata[, !(names(heatdata) %in% c("ext_gene", "log2FoldChange", "padj", "pvalue"))]
#heat2data <- heat2data[,-c(1:8)]

heat2data <- heat2data[, !colnames(heat2data) %in% "z12"]


heatdata1 <- arrange(heat2data,by = padj)
heatdata1 <- heatdata1[1:500,]
#rownames(heatdata1) <- heatdata1$ext_gene



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
                   cluster_cols = T)


pheatmap::pheatmap(heatdata3_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "centroid",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T)


pdf(file = "D:/bulk_control_AS/bulk26/figure/bulk26_heatmap_ward2.pdf")
pheatmap::pheatmap(heatdata3_log,
                   lustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2",
                   show_rownames = F, 
                   show_colnames = T,
                   cluster_cols = T)
dev.off()


#ggsave("bulk26_heatmap_ward2.tiff",device = "tiff",dpi = 300,path = "D:/bulk_control_AS/bulk26/figure")





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

