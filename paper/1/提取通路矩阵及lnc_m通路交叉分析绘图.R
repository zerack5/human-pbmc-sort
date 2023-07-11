



library(KEGGREST)
pathway_list <- keggList("pathway", "hsa")

library(clusterProfiler)

kegg_pathways <- getKEGGPathway()

library(jsonlite)
jdata <- fromJSON("D:/database/GSEA/C2KEGG/c2.cp.kegg.v2023.1.Hs.json")

# 初始化一个空的数据框
result_df <- data.frame(matrix(ncol = 2, nrow = 0))

# 循环遍历largelist的每个元素，提取两个小元素并添加到数据框中
for (i in 1:length(jdata)) {
  sublist <- jdata[[i]]
  row <- c(names(jdata[i]), sublist$exactSource)
  result_df <- rbind(result_df, row)
}

# 设置数据框的列名
colnames(result_df) <- c("pathwayname", "ID")

write.csv(result_df,file = "D:/database/GSEA/C2KEGG/GSEAKEGGHSALIST.csv")


ncCD4c_up <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_cd4_sig_up.rdata")
ncCD4c_down <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_cd4_sig_down.rdata")
ncCD8c_up <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_cd8_sig_up.rdata")
ncCD8c_down <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_cd8_sig_down.rdata")
ncMoc_up <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_mo_sig_up.rdata")
ncMoc_down <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_mo_sig_down.rdata")
ncNKc_up <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_nk_sig_up.rdata")
ncNKc_down <- readRDS("D:/data_anlysis/asncRNA/nckeggresult/2500geneinput/de_nk_sig_down.rdata")

m_CD4_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD4.xlsx",sheet = "KEGG_up")
m_CD4_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD4.xlsx",sheet = "KEGG_down")
m_CD8_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD8.xlsx",sheet = "KEGG_up")
m_CD8_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD8.xlsx",sheet = "KEGG_down")
m_Mo_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_Mo.xlsx",sheet = "KEGG_up")
m_Mo_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_Mo.xlsx",sheet = "KEGG_down")
m_NK_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_Nk.xlsx",sheet = "KEGG_up")
m_NK_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_NK.xlsx",sheet = "KEGG_down")

library(VennDiagram)

venn.diagram(
  x = list(ncCD4c_down@ID, m_CD4_up$ID[1:24]),
  category.names = c("CD4up_mRNA", "CD4up_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD4_lncRNA_mRNAup.png',
  output=TRUE
)
CD4u <- intersect(ncCD4c_down@ID, m_CD4_up$ID[1:24])

venn.diagram(
  x = list(ncCD4c_up@ID, m_CD4_down$ID[1:24]),
  category.names = c("CD4down_mRNA", "CD4down_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD4_lncRNA_mRNAdown.png',
  output=TRUE
)

CD4d <- intersect(ncCD4c_up@ID, m_CD4_down$ID[1:24])

###########################CD8
venn.diagram(
  x = list(ncCD8c_down@ID, m_CD8_up$ID[1:24]),
  category.names = c("CD8up_mRNA", "CD8up_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD8_lncRNA_mRNAup.png',
  output=TRUE
)
CD8u <- intersect(ncCD8c_down@ID, m_CD8_up$ID[1:24])

venn.diagram(
  x = list(ncCD8c_up@ID, m_CD8_down$ID[1:24]),
  category.names = c("CD8down_mRNA", "CD8down_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD8_lncRNA_mRNAdown.png',
  output=TRUE
)

CD8d <- intersect(ncCD8c_up@ID, m_CD8_down$ID[1:24])


###########################Mo
venn.diagram(
  x = list(ncMoc_down@ID, m_Mo_up$ID[1:24]),
  category.names = c("Moup_mRNA", "Moup_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_Mo_lncRNA_mRNAup.png',
  output=TRUE
)
MOu <- intersect(ncMoc_down@ID, m_Mo_up$ID[1:24])

venn.diagram(
  x = list(ncMoc_up@ID, m_Mo_down$ID[1:24]),
  category.names = c("Modown_mRNA", "Modown_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_Mo_lncRNA_mRNAdown.png',
  output=TRUE
)

Mod <- intersect(ncMoc_up@ID, m_Mo_down$ID[1:24])



###########################Nk
venn.diagram(
  x = list(ncNKc_down@ID, m_NK_up$ID[1:24]),
  category.names = c("NKup_mRNA", "NKup_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_NK_lncRNA_mRNAup.png',
  output=TRUE
)
NKu <- intersect(ncNKc_down@ID, m_NK_up$ID[1:24])

venn.diagram(
  x = list(ncNKc_up@ID, m_NK_down$ID[1:24]),
  category.names = c("NKdown_mRNA", "NKdown_lncRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_NK_lncRNA_mRNAdown.png',
  output=TRUE
)

NKd <- intersect(ncNKc_up@ID, m_NK_down$ID[1:24])


########因为lnc通路是

lnc_CD4_up <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/ResultCD4up.rdata")
lnc_CD4_down <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/ResultCD4down.rdata")
lnc_CD8_up <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/ResultCD8up.rdata")
lnc_CD8_down <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/ResultCD8down.rdata")
lnc_mo_up <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/Resultmoup.rdata")
lnc_mo_down <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/Resultmodown.rdata")
lnc_nk_up <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/Resultnkup.rdata")
lnc_nk_down <- readRDS("D:/data_anlysis/asncRNA/LncPathresult/Resultnkdown.rdata")

m_CD4_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD4.xlsx",sheet = "KEGG_up")
m_CD4_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD4.xlsx",sheet = "KEGG_down")
m_CD8_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD8.xlsx",sheet = "KEGG_up")
m_CD8_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_CD8.xlsx",sheet = "KEGG_down")
m_Mo_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_Mo.xlsx",sheet = "KEGG_up")
m_Mo_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_Mo.xlsx",sheet = "KEGG_down")
m_NK_up<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_Nk.xlsx",sheet = "KEGG_up")
m_NK_down<- readxl::read_xlsx("D:/as分选课题/分选课题/KEGG_NK.xlsx",sheet = "KEGG_down")


# 创建一个空的列表用于存储筛选结果

largelist <- lnc_CD4_up
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_CD4_u <- dataframe1

largelist <- lnc_CD4_down
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_CD4_d <- dataframe1


largelist <- lnc_CD8_up
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_CD8_u <- dataframe1

largelist <- lnc_CD8_down
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_CD8_d <- dataframe1


largelist <- lnc_mo_up
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_mo_u <- dataframe1

largelist <- lnc_mo_down
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_mo_d <- dataframe1


largelist <- lnc_nk_up
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_nk_u <- dataframe1


largelist <- lnc_nk_down
filtered_list <- character()
# 使用循环遍历largelist中的每个list
for (i in seq_along(largelist)) {
  # 对每个list进行筛选操作，例如筛选元素大于某个阈值
  filtered <- largelist[[i]][largelist[[i]]$P_Val < 0.01]$gs.name
  # 将筛选结果添加到filtered_list中
  # 检查筛选结果的长度
  if (length(filtered) > 0) {
    # 将筛选结果添加到filtered_list中
    filtered_list[i] <- filtered
  }
  
}
# 删除filtered_list中的NULL值
filtered_list <- na.omit(filtered_list)
dataframe1 <- data.frame(pathwayname = filtered_list)
dataframe1 <- merge(dataframe1,result_df, by ="pathwayname")
lnc_nk_d <- dataframe1









colnames(result_df)[2] <- "ID"
library(VennDiagram)

venn.diagram(
  x = list(lnc_CD4_u$ID,m_CD4_up$ID[1:24]),
  category.names = c( "CD4up_lncRNA","CD4up_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD4_lncRNA_mRNAup.png',
  output=TRUE
)
std <- intersect(lnc_CD4_u$ID, m_CD4_up$ID[1:24])
std <- data.frame(ID = std)
CD4u <- merge(std,result_df, by ="ID")


venn.diagram(
  x = list(lnc_CD4_d$ID,m_CD4_down$ID[1:24]),
  category.names = c( "CD4down_lncRNA","CD4down_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD4_lncRNA_mRNAdown.png',
  output=TRUE
)

std <- intersect(lnc_CD4_d$ID, m_CD4_down$ID[1:24])
std <- data.frame(ID = std)
CD4d <- merge(std,result_df, by ="ID")
###########################CD8
venn.diagram(
  x = list(lnc_CD8_u$ID,m_CD8_up$ID[1:24]),
  category.names = c( "CD8up_lncRNA","CD8up_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD8_lncRNA_mRNAup.png',
  output=TRUE
)
std <- intersect(lnc_CD8_u$ID, m_CD8_up$ID[1:24])
std <- data.frame(ID = std)
CD8u <- merge(std,result_df, by ="ID")

venn.diagram(
  x = list(lnc_CD8_d$ID, m_CD8_down$ID[1:24]),
  category.names = c( "CD8down_lncRNA","CD8down_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_CD8_lncRNA_mRNAdown.png',
  output=TRUE
)

std <- intersect(lnc_CD8_d$ID, m_CD8_down$ID[1:24])
std <- data.frame(ID = std)
CD8d <- merge(std,result_df, by ="ID")

###########################Mo
venn.diagram(
  x = list(lnc_mo_u$ID,m_Mo_up$ID[1:24]),
  category.names = c("Moup_lncRNA","Moup_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_Mo_lncRNA_mRNAup.png',
  output=TRUE
)
std <- intersect(lnc_mo_u$ID, m_Mo_up$ID[1:24])
std <- data.frame(ID = std)
mou <- merge(std,result_df, by ="ID")


venn.diagram(
  x = list(lnc_mo_d$ID,m_Mo_down$ID[1:24]),
  category.names = c("Modown_lncRNA","Modown_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_Mo_lncRNA_mRNAdown.png',
  output=TRUE
)

std <- intersect(lnc_mo_d$ID, m_Mo_down$ID[1:24])
std <- data.frame(ID = std)
mod <- merge(std,result_df, by ="ID")


###########################Nk
venn.diagram(
  x = list(lnc_nk_u$ID,m_NK_up$ID[1:24]),
  category.names = c("NKup_lncRNA","NKup_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_NK_lncRNA_mRNAup.png',
  output=TRUE
)
std <- intersect(lnc_nk_u$ID, m_NK_up$ID[1:24])
std <- data.frame(ID = std)
nku <- merge(std,result_df, by ="ID")

venn.diagram(
  x = list(lnc_nk_d$ID,m_NK_down$ID[1:24]),
  category.names = c("NKdown_lncRNA","NKdown_mRNA"),
  fill = c( "purple","orange"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c( "purple","orange"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = 'D:/as分选课题/分选课题/figure/韦恩图/sorted_NK_lncRNA_mRNAdown.png',
  output=TRUE
)

std <- intersect(lnc_nk_d$ID, m_NK_down$ID[1:24])
std <- data.frame(ID = std)
nkd <- merge(std,result_df, by ="ID")
















names(jdata)
jdata[[3]]
