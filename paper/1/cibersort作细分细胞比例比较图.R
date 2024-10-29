library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
heatdata<-read.csv("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/expression_matrix/Moasssort.csv", header = TRUE)
colnames(heatdata)
heatdata <- heatdata[,c("ext_gene","log2FoldChange","padj","pvalue" )]
#heatdata <- heatdata[which(abs(heatdata$log2FoldChange)>1),]

library(ggplot2)

#将矩阵根据列名中是否含有c一分为2
# 指定要查找的字符
result1 <- read.csv("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/expression_matrix/Mocontrolsort.csv")
result2 <- read.csv("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/expression_matrix/Moasssort.csv")

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

#####处理一下读取的矩阵元素
rownames(result1) <-result1$X
result1$X <- NULL

rownames(result2) <-result2$X
result2$X <- NULL

resultall <- rbind(result1,result2)

##########486版细胞比例


res <- data.frame(resultall[,1:(ncol(resultall) - 3)])%>%
  mutate(group = c(rep('control',nrow(result1)),rep('as',nrow(result2))))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:(ncol(resultall)-2)],
               names_to = "cell.type",
               values_to = 'value')
head(res)
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


###################################################################
# 筛选 p 值显著的结果
#filtered_res <- filter(res, pvalue <= 0.05)
filtered_res <- res
###########################################################


# 根据p值添加显著性标记
filtered_res <- filtered_res %>% 
  mutate(significance = case_when(
    pvalue <= 0.001 ~ "***",
    pvalue <= 0.01 ~ "**",
    pvalue <= 0.05 ~ "*",
    #TRUE ~ NA_character_ 
    #TRUE ~ "ns"
    is.na(pvalue) | pvalue > 0.05 ~ "ns"  # 确保所有不显著的 p 值标记为 "ns"
  ))

# 绘制带显著性标记的箱线图
ggplot(filtered_res, aes(cell.type, value, fill = group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion", size = 10) +
  theme(legend.position = "top",
        legend.text = element_text(size = 15),  # 修改图例文本的大小
        legend.title = element_text(size = 18),  # 修改图例标题的大小
        legend.key.size = unit(1.5, "cm"),       # 修改图例键（颜色块）的大小
        axis.text.x = element_text(angle = 80, vjust = 0.5, size = 15, face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic", size = 15, colour = 'black'),
        axis.title.x = element_text(size = 20),  # 修改横坐标标题大小
        axis.title.y = element_text(size = 20),  # 修改纵坐标标题大小
        panel.grid.major = element_blank(),  # 去除主网格线
        panel.grid.minor = element_blank()   # 去除次网格线
  ) +
  scale_fill_manual(name = "Group",
                    values = c("control" = "steelblue", "as" = "darkred"),
                    labels = c("control" = "Control", "as" = "AS")  # 自定义图例标签
  ) +
  geom_text(data = filtered_res, aes(x = cell.type, y = max(filtered_res$value) + 0.01, label = significance), 
            size = 6, color = "black", vjust = -0.1)
ggsave("PBMC26_lnc_486_cibersort_mo.tiff",device = "tiff",width =18,height =9,dpi = 300,path = "E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/newfigure")




#将矩阵根据列名中是否含有c一分为2
# 指定要查找的字符
result1 <- read.csv("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/expression_matrix/CD4controlsort.csv")
result2 <- read.csv("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/expression_matrix/CD4asssort.csv")

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

#####处理一下读取的矩阵元素
rownames(result1) <-result1$X
result1$X <- NULL

rownames(result2) <-result2$X
result2$X <- NULL

resultall <- rbind(result1,result2)

##########486版细胞比例


res <- data.frame(resultall[,1:(ncol(resultall) - 3)])%>%
  mutate(group = c(rep('control',nrow(result1)),rep('as',nrow(result2))))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:(ncol(resultall)-2)],
               names_to = "cell.type",
               values_to = 'value')
head(res)
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


###################################################################
# 筛选 p 值显著的结果
#filtered_res <- filter(res, pvalue <= 0.05)
filtered_res <- res
###########################################################


# 根据p值添加显著性标记
filtered_res <- filtered_res %>% 
  mutate(significance = case_when(
    pvalue <= 0.001 ~ "***",
    pvalue <= 0.01 ~ "**",
    pvalue <= 0.05 ~ "*",
    #TRUE ~ NA_character_ 
    #TRUE ~ "ns"
    is.na(pvalue) | pvalue > 0.05 ~ "ns"  # 确保所有不显著的 p 值标记为 "ns"
  ))

# 绘制带显著性标记的箱线图
ggplot(filtered_res, aes(cell.type, value, fill = group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion", size = 10) +
  theme(legend.position = "top",
        legend.text = element_text(size = 15),  # 修改图例文本的大小
        legend.title = element_text(size = 18),  # 修改图例标题的大小
        legend.key.size = unit(1.5, "cm"),       # 修改图例键（颜色块）的大小
        axis.text.x = element_text(angle = 80, vjust = 0.5, size = 15, face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic", size = 15, colour = 'black'),
        axis.title.x = element_text(size = 20),  # 修改横坐标标题大小
        axis.title.y = element_text(size = 20),  # 修改纵坐标标题大小
        panel.grid.major = element_blank(),  # 去除主网格线
        panel.grid.minor = element_blank()   # 去除次网格线
  ) +
  scale_fill_manual(name = "Group",
                    values = c("control" = "steelblue", "as" = "darkred"),
                    labels = c("control" = "Control", "as" = "AS")  # 自定义图例标签
  ) +
  geom_text(data = filtered_res, aes(x = cell.type, y = max(filtered_res$value) + 0.01, label = significance), 
            size = 6, color = "black", vjust = -0.1)
ggsave("PBMC26_lnc_486_cibersort_CD4.tiff",device = "tiff",width =18,height =9,dpi = 300,path = "E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/newfigure")






#将矩阵根据列名中是否含有c一分为2
# 指定要查找的字符
result1 <- read.csv("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/expression_matrix/CD8controlsort.csv")
result2 <- read.csv("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/expression_matrix/CD8asssort.csv")

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

#####处理一下读取的矩阵元素
rownames(result1) <-result1$X
result1$X <- NULL

rownames(result2) <-result2$X
result2$X <- NULL

resultall <- rbind(result1,result2)

##########486版细胞比例


res <- data.frame(resultall[,1:(ncol(resultall) - 3)])%>%
  mutate(group = c(rep('control',nrow(result1)),rep('as',nrow(result2))))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:(ncol(resultall)-2)],
               names_to = "cell.type",
               values_to = 'value')
head(res)
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


###################################################################
# 筛选 p 值显著的结果
#filtered_res <- filter(res, pvalue <= 0.05)
filtered_res <- res
###########################################################


# 根据p值添加显著性标记
filtered_res <- filtered_res %>% 
  mutate(significance = case_when(
    pvalue <= 0.001 ~ "***",
    pvalue <= 0.01 ~ "**",
    pvalue <= 0.05 ~ "*",
    #TRUE ~ NA_character_ 
    #TRUE ~ "ns"
    is.na(pvalue) | pvalue > 0.05 ~ "ns"  # 确保所有不显著的 p 值标记为 "ns"
  ))

# 绘制带显著性标记的箱线图
ggplot(filtered_res, aes(cell.type, value, fill = group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion", size = 10) +
  theme(legend.position = "top",
        legend.text = element_text(size = 15),  # 修改图例文本的大小
        legend.title = element_text(size = 18),  # 修改图例标题的大小
        legend.key.size = unit(1.5, "cm"),       # 修改图例键（颜色块）的大小
        axis.text.x = element_text(angle = 80, vjust = 0.5, size = 15, face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic", size = 15, colour = 'black'),
        axis.title.x = element_text(size = 20),  # 修改横坐标标题大小
        axis.title.y = element_text(size = 20),  # 修改纵坐标标题大小
        panel.grid.major = element_blank(),  # 去除主网格线
        panel.grid.minor = element_blank()   # 去除次网格线
  ) +
  scale_fill_manual(name = "Group",
                    values = c("control" = "steelblue", "as" = "darkred"),
                    labels = c("control" = "Control", "as" = "AS")  # 自定义图例标签
  ) +
  geom_text(data = filtered_res, aes(x = cell.type, y = max(filtered_res$value) + 0.01, label = significance), 
            size = 6, color = "black", vjust = -0.1)
ggsave("PBMC26_lnc_486_cibersort_CD8.tiff",device = "tiff",width =18,height =9,dpi = 300,path = "E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/newfigure")




