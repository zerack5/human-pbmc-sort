######韦恩图
setwd("D:/data_anlysis/asncRNA/LncPathresult/FDRcut0.01")
library(VennDiagram)

CD4up <- read.table("unique_CD4up03.txt")
CD4up <- CD4up$V1

CD8up <- read.table("unique_CD8up03.txt")
CD8up <- CD8up$V1

moup <- read.table("unique_moup03.txt")
moup <- moup$V1

nkup <- read.table("unique_nkup03.txt")
nkup <- nkup$V1


venn.diagram(
  x = list(CD4up, CD8up, moup, nkup),
  category.names = c("CD4up", "CD8up", "moup","nkup"),
  fill = c("skyblue", "mediumorchid", "lightgreen","yellow"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c("skyblue", "mediumorchid", "lightgreen","yellow"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = '#1_sorted_lncRNA_up.png',
  output=TRUE
)

CD4down <- read.table("unique_CD4down03.txt")
CD4down <- CD4down$V1

CD8down <- read.table("unique_CD8down03.txt")
CD8down <- CD8down$V1

modown <- read.table("unique_modown03.txt")
modown <- modown$V1

nkdown <- read.table("unique_nkdown03.txt")
nkdown <- nkdown$V1


venn.diagram(
  x = list(CD4down, CD8down, modown, nkdown),
  category.names = c("CD4down", "CD8down", "modown","nkdown"),
  fill = c("skyblue", "mediumorchid", "lightgreen","yellow"),
  alpha = 0.5,
  label.col = "black",
  cex = 0.8,
  cat.col = c("skyblue", "mediumorchid", "lightgreen","yellow"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  filename = '#1_sorted_lncRNA_down.png',
  output=TRUE
)




# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)