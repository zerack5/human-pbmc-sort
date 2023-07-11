#########asnc_bulkdata

path0 <- "D:/data_anlysis/asncRNA/cibersort"


####构建细分亚群的标记矩阵


setwd(path0)
hc486 <- read.table("D:/data_anlysis/public_data_of_the_cell_paper/cibersortx/486hcreferrence_sample_file.txt", header = T, row.names = 1)

row.names(hc486)<- hc486$Gene

#####CD4T细胞亚群
# 指定要筛选的特定字符串
target_strings <- c("Naive_CD4", "Mem_CD4","Th1","Th2","Th17","Tfh","Fr._I_nTreg","Fr._II_eTreg","Fr._III_T")
# 提取包含特定字符串的列名
matching_cols <- grepl(paste(target_strings, collapse = "|"), names(hc486))
# 提取含有特定字符串的列形成新的数据框
CD4hc486 <- hc486[,matching_cols]

merged_df <- CD4hc486
stmatrix <- matrix(nrow = length(target_strings),ncol = ncol( merged_df)-1)
colnames(stmatrix) <- colnames(merged_df)[-1]
rownames(stmatrix) <- target_strings
result <- stmatrix
# 遍历df中的每一个元素
for (i in 1:nrow(stmatrix)) {
  for (j in 1:ncol(stmatrix)) {
    # 比较行名与部分列名,使用.号或者下划线之前的列名比较
    if (rownames(stmatrix)[i] == sub("_[a-z0-9A-Z]+$","",colnames(stmatrix))[j]) {
      result[i, j] <- 1
    } else {
      result[i, j] <- 2
    }
  }
}
t_merged_df <- merged_df
colnames(t_merged_df) <- sub("^(.*)_.*$", "\\1", colnames(t_merged_df))
write.table(t_merged_df, file = "D:/data_anlysis/asncRNA/cibersort/sample_file/CD4_486hcreferrence_sample_file.txt", sep = "\t",quote = F,col.names = T,row.names = T)
write.table(result, file = "D:/data_anlysis/asncRNA/cibersort/phenotype/CD4_486hcphenotype_classes_file.txt", sep = "\t",col.names = F,quote = F)






######CD8T细胞亚群
target_strings <- c("Naive_CD8", "CM_CD8","EM_CD8","TEMRA_CD8")
# 提取包含特定字符串的列名
matching_cols <- grepl(paste(target_strings, collapse = "|"), names(hc486))
# 提取含有特定字符串的列形成新的数据框
CD8hc486 <- hc486[,matching_cols]

merged_df <- CD8hc486
stmatrix <- matrix(nrow = length(target_strings),ncol = ncol( merged_df)-1)
colnames(stmatrix) <- colnames(merged_df)[-1]
rownames(stmatrix) <- target_strings
result <- stmatrix
# 遍历df中的每一个元素
for (i in 1:nrow(stmatrix)) {
  for (j in 1:ncol(stmatrix)) {
    # 比较行名与部分列名,使用.号或者下划线之前的列名比较
    if (rownames(stmatrix)[i] == sub("_[a-z0-9A-Z]+$","",colnames(stmatrix))[j]) {
      result[i, j] <- 1
    } else {
      result[i, j] <- 2
    }
  }
}
t_merged_df <- merged_df
colnames(t_merged_df) <- sub("^(.*)_.*$", "\\1", colnames(t_merged_df))
write.table(t_merged_df, file = "D:/data_anlysis/asncRNA/cibersort/sample_file/CD8_486hcreferrence_sample_file.txt", sep = "\t",quote = F,col.names = T,row.names = T)
write.table(result, file = "D:/data_anlysis/asncRNA/cibersort/phenotype/CD8_486hcphenotype_classes_file.txt", sep = "\t",col.names = F,quote = F)







#########Monocyte
target_strings <- c("CL_Mono", "CD16p_Mono","Int_Mono","NC_Mono")
# 提取包含特定字符串的列名
matching_cols <- grepl(paste(target_strings, collapse = "|"), names(hc486))
# 提取含有特定字符串的列形成新的数据框
Mohc486 <- hc486[,matching_cols]

merged_df <- Mohc486
stmatrix <- matrix(nrow = length(target_strings),ncol = ncol( merged_df)-1)
colnames(stmatrix) <- colnames(merged_df)[-1]
rownames(stmatrix) <- target_strings
result <- stmatrix
# 遍历df中的每一个元素
for (i in 1:nrow(stmatrix)) {
  for (j in 1:ncol(stmatrix)) {
    # 比较行名与部分列名,使用.号或者下划线之前的列名比较
    if (rownames(stmatrix)[i] == sub("_[a-z0-9A-Z]+$","",colnames(stmatrix))[j]) {
      result[i, j] <- 1
    } else {
      result[i, j] <- 2
    }
  }
}
t_merged_df <- merged_df
colnames(t_merged_df) <- sub("^(.*)_.*$", "\\1", colnames(t_merged_df))
write.table(t_merged_df, file = "D:/data_anlysis/asncRNA/cibersort/sample_file/Mo_486hcreferrence_sample_file.txt", sep = "\t",quote = F,col.names = T,row.names = T)
write.table(result, file = "D:/data_anlysis/asncRNA/cibersort/phenotype/Mo_486hcphenotype_classes_file.txt", sep = "\t",col.names = F,quote = F)







#####CD4T细胞亚群






######CD8T细胞亚群




#####Mo细胞亚群




#####对分选后的细胞进行构图












