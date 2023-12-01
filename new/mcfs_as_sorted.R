rm(list = ls())
options(java.parameters = "-Xmx4g")
library(rmcfs)
library(ggplot2)


setwd("C:/Users/a/Desktop/表达矩阵")

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



write.csv(heatCD4data,file="D:/bulk_control_AS/svm/cd4_count.csv")
write.csv(heatCD8data,file="D:/bulk_control_AS/svm/cd8_count.csv")
write.csv(heatmodata,file="D:/bulk_control_AS/svm/mo_count.csv")
write.csv(heatnkdata,file="D:/bulk_control_AS/svm/nk_count.csv")







# 定义一个函数来根据列名的模式添加新的行
add_control_or_AS <- function(data) {
  # 获取数据框的列名
  col_names <- colnames(data)
  
  # 创建新的行
  new_row <- numeric(length(col_names))
  
  # 遍历列名
  for (i in 1:length(col_names)) {
    if (grepl("^C\\d+", col_names[i])) {  # 列名匹配C紧接着任意数字
      new_row[i] <- "Control"
    } else {
      new_row[i] <- "AS"
    }
  }
  
  # 将新行添加到数据框末尾
  data[nrow(data) + 1, ] <- new_row
  
  # 更改新行的行名为"label"
  rownames(data)[nrow(data)] <- "label" 
  
  
  return(data)
}

# 调用函数来处理数据框
#new_df <- add_control_or_AS(df)



CD4data <- read.csv("D:/bulk_control_AS/svm/cd4_count.csv",row.names = 2)
CD4data$X <- NULL
new_CD4data <- add_control_or_AS(CD4data)
new_CD4data <- t(new_CD4data)
new_CD4data <- as.data.frame(new_CD4data)


CD8data <- read.csv("D:/bulk_control_AS/svm/cd8_count.csv",row.names = 2)
CD8data$X <- NULL
new_CD8data <- add_control_or_AS(CD8data)
new_CD8data <- t(new_CD8data)
new_CD8data <- as.data.frame(new_CD8data)


modata <- read.csv("D:/bulk_control_AS/svm/mo_count.csv",row.names = 2)
modata$X <- NULL
new_modata <- add_control_or_AS(modata)
new_modata <- t(new_modata)
new_modata <- as.data.frame(new_modata)



nkdata <- read.csv("D:/bulk_control_AS/svm/nk_count.csv",row.names = 2)
nkdata$X <- NULL
new_nkdata <- add_control_or_AS(nkdata)
new_nkdata <- t(new_nkdata)
new_nkdata <- as.data.frame(new_nkdata)






#z_smote_CRC_2 <- read.table("smote_CRC_no_zscore.txt", header = T, sep = "\t")
# Fix data types and data values - replace characters such as "," " " "/" etc. 
# from values and column names and fix data types
# This function may help if mcfs has any problems with input data

#CRC_mcfs <- fix.data(z_smote_CRC_2)
# Run MCFS-ID procedure on default parameters. 
# For larger real data (thousands of features) default 'auto' settings are the best.
# This example may take 10-20 minutes but this one is a real dataset with 4026 features.
# Set up more threads according to your CPU cores number.



new_testdata <- fix.data(new_CD4data)
setwd("D:/bulk_control_AS/svm")
result_CRC <- mcfs(label~., 
                   new_testdata, 
                   featureFreq = 100, 
                   cutoffMethod = c("permutations", "criticalAngle", "kmeans", "mean", "contrast"),
                   cutoffPermutations = 20, 
                   finalCV = TRUE,
                   seed = 2022, 
                   threadsNumber = 8)
save(result_CRC, file = "D:/bulk_control_AS/svm/cd4/all_result_CD4_1.RData")
load("D:/bulk_control_AS/svm/cd4/all_result_CD4_1.RData")

# Print basic information about mcfs result.
all_cutoff <- result_CRC$cutoff
pdf("cd4/cd4_mcfs_size.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = log10(size))) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

all_minRI <- result_CRC$cutoff
pdf("cd4/cd4_mcfs_minRI.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = minRI)) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

# 结合上面两张图之后，选择permutations，不需要再用permutation再做一次，验证过了

# Plot cv classification result obtained on top features.
# In the middle of x axis red label denotes cutoff_value.
pdf("cd4/cd4_cv_acc.pdf",width = 5, height = 5)
plot(result_CRC, type = "cv", cv_measure = "acc", cex = 0.8)
dev.off()

# Plot & print out confusion matrix. This matrix is the result of 
# all classifications performed by all decision trees on all s*t datasets.
pdf("cd4/cd4_mcfs_cmatrix.pdf", width = 10)
plot(result_CRC, type = "cmatrix")
dev.off()

mcfs_RI <- result_CRC$RI;dim(mcfs_RI)
mcfs_name <- mcfs_RI$attribute
write.table(mcfs_RI, "cd4/all_mcfs_cd4_rank.txt", sep = "\t", quote = F, row.names = F)
write.table(new_testdata[, c(mcfs_name, "label")], "cd4/all_mcfs_rank_expr.txt", 
            sep = "\t", quote = F, row.names = F)



############################Cd8
#####这一步fix过了以后改变了列名

new_testdata <- fix.data(new_CD8data)
setwd("D:/bulk_control_AS/svm")
result_CRC <- mcfs(label~., 
                   new_testdata, 
                   featureFreq = 100, 
                   cutoffMethod = c("permutations", "criticalAngle", "kmeans", "mean", "contrast"),
                   cutoffPermutations = 20, 
                   finalCV = TRUE,
                   seed = 2022, 
                   threadsNumber = 8)
save(result_CRC, file = "D:/bulk_control_AS/svm/cd8/all_result_CD8_1.RData")
load("D:/bulk_control_AS/svm/cd8/all_result_CD8_1.RData")

# Print basic information about mcfs result.
all_cutoff <- result_CRC$cutoff
pdf("cd8/cd8_mcfs_size.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = log10(size))) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

all_minRI <- result_CRC$cutoff
pdf("cd8/cd8_mcfs_minRI.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = minRI)) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

# 结合上面两张图之后，选择permutations，不需要再用permutation再做一次，验证过了

# Plot cv classification result obtained on top features.
# In the middle of x axis red label denotes cutoff_value.
pdf("cd8/cd8_cv_acc.pdf",width = 5, height = 5)
plot(result_CRC, type = "cv", cv_measure = "acc", cex = 0.8)
dev.off()

# Plot & print out confusion matrix. This matrix is the result of 
# all classifications performed by all decision trees on all s*t datasets.
pdf("cd8/cd8_mcfs_cmatrix.pdf", width = 10)
plot(result_CRC, type = "cmatrix")
dev.off()

mcfs_RI <- result_CRC$RI;dim(mcfs_RI)
mcfs_name <- mcfs_RI$attribute
write.table(mcfs_RI, "cd8/all_mcfs_cd8_rank.txt", sep = "\t", quote = F, row.names = F)
write.table(new_testdata[, c(mcfs_name, "label")], "cd8/cd8_mcfs_rank_expr.txt", 
            sep = "\t", quote = F, row.names = F)
#####这一步fix过了以后改变了列名


#####################################################mo



new_testdata <- fix.data(new_modata)
setwd("D:/bulk_control_AS/svm")
result_CRC <- mcfs(label~., 
                   new_testdata, 
                   featureFreq = 100, 
                   cutoffMethod = c("permutations", "criticalAngle", "kmeans", "mean", "contrast"),
                   cutoffPermutations = 20, 
                   finalCV = TRUE,
                   seed = 2022, 
                   threadsNumber = 8)
save(result_CRC, file = "D:/bulk_control_AS/svm/mo/all_result_mo_1.RData")
load("D:/bulk_control_AS/svm/mo/all_result_mo_1.RData")

# Print basic information about mcfs result.
all_cutoff <- result_CRC$cutoff
pdf("mo/mo_mcfs_size.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = log10(size))) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

all_minRI <- result_CRC$cutoff
pdf("mo/mo_mcfs_minRI.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = minRI)) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

# 结合上面两张图之后，选择permutations，不需要再用permutation再做一次，验证过了

# Plot cv classification result obtained on top features.
# In the middle of x axis red label denotes cutoff_value.
pdf("mo/mo_cv_acc.pdf",width = 5, height = 5)
plot(result_CRC, type = "cv", cv_measure = "acc", cex = 0.8)
dev.off()

# Plot & print out confusion matrix. This matrix is the result of 
# all classifications performed by all decision trees on all s*t datasets.
pdf("mo/mo_mcfs_cmatrix.pdf", width = 10)
plot(result_CRC, type = "cmatrix")
dev.off()

mcfs_RI <- result_CRC$RI;dim(mcfs_RI)
mcfs_name <- mcfs_RI$attribute
write.table(mcfs_RI, "mo/all_mcfs_mo_rank.txt", sep = "\t", quote = F, row.names = F)
write.table(new_testdata[, c(mcfs_name, "label")], "mo/all_mcfs_rank_expr.txt", 
            sep = "\t", quote = F, row.names = F)



###################################################################nk




new_testdata <- fix.data(new_nkdata)
setwd("D:/bulk_control_AS/svm")
result_CRC <- mcfs(label~., 
                   new_testdata, 
                   featureFreq = 100, 
                   cutoffMethod = c("permutations", "criticalAngle", "kmeans", "mean", "contrast"),
                   cutoffPermutations = 20, 
                   finalCV = TRUE,
                   seed = 2022, 
                   threadsNumber = 8)
save(result_CRC, file = "D:/bulk_control_AS/svm/nk/all_result_nk_1.RData")
load("D:/bulk_control_AS/svm/nk/all_result_nk_1.RData")

# Print basic information about mcfs result.
all_cutoff <- result_CRC$cutoff
pdf("nk/nk_mcfs_size.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = log10(size))) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

all_minRI <- result_CRC$cutoff
pdf("nk/nk_mcfs_minRI.pdf", width = 6, height = 6)
ggplot(data = all_cutoff, aes(x = method, y = minRI)) +
  geom_bar(stat = 'identity',
           fill = c("#ee96dc", "#92cee2", "#ffd977", "#ff9747"),
           color = "black") +
  theme_test() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0))
dev.off()

# 结合上面两张图之后，选择permutations，不需要再用permutation再做一次，验证过了

# Plot cv classification result obtained on top features.
# In the middle of x axis red label denotes cutoff_value.
pdf("nk/nk_cv_acc.pdf",width = 5, height = 5)
plot(result_CRC, type = "cv", cv_measure = "acc", cex = 0.8)
dev.off()

# Plot & print out confusion matrix. This matrix is the result of 
# all classifications performed by all decision trees on all s*t datasets.
pdf("nk/nk_mcfs_cmatrix.pdf", width = 10)
plot(result_CRC, type = "cmatrix")
dev.off()

mcfs_RI <- result_CRC$RI;dim(mcfs_RI)
mcfs_name <- mcfs_RI$attribute
write.table(mcfs_RI, "nk/all_mcfs_nk_rank.txt", sep = "\t", quote = F, row.names = F)
write.table(new_testdata[, c(mcfs_name, "label")], "nk/all_mcfs_rank_expr.txt", 
            sep = "\t", quote = F, row.names = F)






