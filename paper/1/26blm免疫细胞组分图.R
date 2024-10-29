


















#将矩阵根据列名中是否含有c一分为2
# 指定要查找的字符#
#这里的命名错了，懒得改了标注一下
result1 <- read.csv("LM22_casecontrol_cibersortresult.csv")
result2 <- read.csv("Cell486_casecontrol_cibersortresult.csv")

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


#cibersort细胞比例分析

######自建细胞矩阵

# 数据准备
res <- data.frame(result1[, 1:27]) %>%
  mutate(group = c(rep('control', 26), rep('as', 26))) %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = colnames(.)[2:28],
               names_to = "cell.type",
               values_to = 'value')
res$group <- factor(res$group, levels = c('control', 'as'))

# 计算 p 值
res$pvalue <- NA
for (i in unique(res$cell.type)) {
  the_dt <- res[which(res$cell.type == i), ]
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
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 80, vjust = 0.5, size = 14, face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic", size = 14, colour = 'black'),
        panel.grid.major = element_blank(),  # 去除主网格线
        panel.grid.minor = element_blank()   # 去除次网格线
        
        ) +
  scale_fill_manual(name= "Group",
                    values = c("control" = "steelblue", "as" = "darkred"),
                    labels = c("control" = "Control", "as" = "AS")  # 自定义图例标签
                    ) +
  geom_text(data = filtered_res, aes(x = cell.type, y = max(filtered_res$value) + 0.01, label = significance), 
            size = 5, color = "black", vjust = -0.5)




ggsave("PBMC26_m_486_cibersort_all.tiff",device = "tiff",width =18,height =9,dpi = 300,path = "E:/同步文件夹/BaiduSyncdisk/data_anlysis/bulk_control_AS/figure96/cibersort/")


#########原版LM22细胞矩阵





# 数据准备
res <- data.frame(result2[,1:22])%>%
  mutate(group = c(rep('control',26),rep('as',26)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value')
head(res,n=42)
res$group <- factor(res$group, levels = c('control', 'as'))

# 计算 p 值
res$pvalue <- NA
for (i in unique(res$cell.type)) {
  the_dt <- res[which(res$cell.type == i), ]
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
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 80, vjust = 0.5, size = 14, face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic", size = 14, colour = 'black'),
        panel.grid.major = element_blank(),  # 去除主网格线
        panel.grid.minor = element_blank()   # 去除次网格线
        
  ) +
  scale_fill_manual(name= "Group",
                    values = c("control" = "steelblue", "as" = "darkred"),
                    labels = c("control" = "Control", "as" = "AS")  # 自定义图例标签
  ) +
  geom_text(data = filtered_res, aes(x = cell.type, y = max(filtered_res$value) + 0.01, label = significance), 
            size = 5, color = "black", vjust = -0.5)




ggsave("PBMC26_m_lm22_cibersort_all.tiff",device = "tiff",width =18,height =9,dpi = 300,path = "E:/同步文件夹/BaiduSyncdisk/data_anlysis/bulk_control_AS/figure96/cibersort/")




result_mo <- read.table("E:/同步文件夹/BaiduSyncdisk/data_anlysis/asncRNA/cibersort/sig_matrix/Mo_hc486_sig_matrix.txt")


























