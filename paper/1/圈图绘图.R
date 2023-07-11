#########圈图
path0 <- "D:/data_anlysis/asncRNA"
path1 <- "D:/data_anlysis/asncRNA/lnc_mrna_coexpression"

cnamecount <- "ncrna_control_count.csv"
asnamecount <- "ncrna_as_sum_count.csv"  
sname1sum <- "t1phase_il17asum.txt"
sname2 <- "t2phase_il17a.txt"
sname2sum <- "t2phase_il17asum.txt"
sname3 <- "t3phase_il17a.txt"
sname3sum <- "t3phase_il17asum.txt"


lncnamec <- "lnc_mrna_coexpression/ncrna_control_genename_count.csv"
lncnameas <- "lnc_mrna_coexpression/ncrna_as_genename_count.csv"

mnamec <- "C:/Users/a/Desktop/Rrun1/exp2control.csv"
mnameas <- "C:/Users/a/Desktop/Rrun1/exp2as.csv"


lnc_deseq_cd4 <- "D:/data_anlysis/asncRNA/deseqresult/deseq2_CD4_sig.csv"
lnc_deseq_cd8 <- "D:/data_anlysis/asncRNA/deseqresult/deseq2_CD8_sig.csv"
lnc_deseq_mo <- "D:/data_anlysis/asncRNA/deseqresult/deseq2_mo_sig.csv"
lnc_deseq_nk <- "D:/data_anlysis/asncRNA/deseqresult/deseq2_nk_sig.csv"

m_deseq_cd4 <- "D:/as分选课题/分选课题/1/CD4_res.csv"
m_deseq_cd8 <- "D:/as分选课题/分选课题/1/CD8_res.csv"
m_deseq_mo <- "D:/as分选课题/分选课题/1/res_Mo.csv"
m_deseq_nk <- "D:/as分选课题/分选课题/1/res_NK.csv"


library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

t7g <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol',
                                 'chromosome_name','start_position','end_position',"strand"), 
                    mart = ensembl)
# 获取属性列表
#attributes <- listAttributes(ensembl)
#martlist <- listMarts()
# 打印属性列表
#print(attributes)


library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

t6g <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',
                          'chromosome_name','start_position','end_position'), 
             mart = ensembl)
saveRDS(t6g, file = "D:/Rdata/biomart/t6g.rdata")

t6g <- readRDS("D:/Rdata/biomart/t6g.rdata")
#save(t7g, file = "D:/Rdata/biomart/t7g.rdata")


lnc_deseqre_cd4 <- read.csv("D:/data_anlysis/asncRNA/deseqresult/deseq2_CD4_sig.csv",header = T)
lnc_deseqre_cd8 <- read.csv("D:/data_anlysis/asncRNA/deseqresult/deseq2_CD8_sig.csv",header = T)
lnc_deseqre_mo <- read.csv("D:/data_anlysis/asncRNA/deseqresult/deseq2_mo_sig.csv",header = T)
lnc_deseqre_nk <- read.csv("D:/data_anlysis/asncRNA/deseqresult/deseq2_nk_sig.csv",header = T)





colnames(lnc_deseqre_cd4)[1] <- "ensembl_gene_id"
colnames(lnc_deseqre_cd8)[1] <- "ensembl_gene_id"
colnames(lnc_deseqre_mo)[1] <- "ensembl_gene_id"
colnames(lnc_deseqre_nk)[1] <- "ensembl_gene_id"
lnc_deseqre_cd4 <- merge(lnc_deseqre_cd4,t6g,by = "ensembl_gene_id")
lnc_deseqre_cd8 <- merge(lnc_deseqre_cd8,t6g,by = "ensembl_gene_id")
lnc_deseqre_mo <- merge(lnc_deseqre_mo,t6g,by = "ensembl_gene_id")
lnc_deseqre_nk <- merge(lnc_deseqre_nk,t6g,by = "ensembl_gene_id")
#将external_gene_name中的空字符串替换成ensemnl_gene_id
lnc_deseqre_cd4$hgnc_symbol[lnc_deseqre_cd4$hgnc_symbol==""] <- lnc_deseqre_cd4$ensembl_gene_id[lnc_deseqre_cd4$hgnc_symbol==""]
lnc_deseqre_cd8$hgnc_symbol[lnc_deseqre_cd8$hgnc_symbol==""] <- lnc_deseqre_cd8$ensembl_gene_id[lnc_deseqre_cd8$hgnc_symbol==""]
lnc_deseqre_mo$hgnc_symbol[lnc_deseqre_mo$hgnc_symbol==""] <- lnc_deseqre_mo$ensembl_gene_id[lnc_deseqre_mo$hgnc_symbol==""]
lnc_deseqre_nk$hgnc_symbol[lnc_deseqre_nk$hgnc_symbol==""] <- lnc_deseqre_nk$ensembl_gene_id[lnc_deseqre_nk$hgnc_symbol==""]
#将显著差异表达的基因名做成一个字符长串
lnc_deseqsg_cd4 <- lnc_deseqre_cd4[abs(lnc_deseqre_cd4$padj)<0.001&abs(lnc_deseqre_cd4$log2FoldChange)>1,]
lnc_deseqsg_cd8 <- lnc_deseqre_cd8[abs(lnc_deseqre_cd8$padj)<0.001&abs(lnc_deseqre_cd8$log2FoldChange)>1,]
lnc_deseqsg_mo <- lnc_deseqre_mo[abs(lnc_deseqre_mo$padj)<0.001&abs(lnc_deseqre_mo$log2FoldChange)>1,]
lnc_deseqsg_nk <- lnc_deseqre_nk[abs(lnc_deseqre_nk$padj)<0.001&abs(lnc_deseqre_nk$log2FoldChange)>1,]

#使用正则表达式和grepl()函数检查列中的每个元素是否为数字或"X"和"Y"字符：
numeric_mask <- grepl("^[0-9XY]+$", lnc_deseqsg_cd4$chromosome_name)
lnc_deseqsg_cd4 <- subset(lnc_deseqsg_cd4, numeric_mask)
lnc_deseqsg_cd4$chromosome_name <- paste0("chr", lnc_deseqsg_cd4$chromosome_name)

numeric_mask <- grepl("^[0-9XY]+$", lnc_deseqsg_cd8$chromosome_name)
lnc_deseqsg_cd8 <- subset(lnc_deseqsg_cd8, numeric_mask)
lnc_deseqsg_cd8$chromosome_name <- paste0("chr", lnc_deseqsg_cd8$chromosome_name)

numeric_mask <- grepl("^[0-9XY]+$", lnc_deseqsg_mo$chromosome_name)
lnc_deseqsg_mo <- subset(lnc_deseqsg_mo, numeric_mask)
lnc_deseqsg_mo$chromosome_name <- paste0("chr", lnc_deseqsg_mo$chromosome_name)

numeric_mask <- grepl("^[0-9XY]+$", lnc_deseqsg_nk$chromosome_name)
lnc_deseqsg_nk <- subset(lnc_deseqsg_nk, numeric_mask)
lnc_deseqsg_nk$chromosome_name <- paste0("chr", lnc_deseqsg_nk$chromosome_name)





m_deseqre_cd4 <- read.csv("C:/Users/a/Desktop/333/all_diffsigcd4.csv",header = T)
m_deseqre_cd8 <- read.csv("C:/Users/a/Desktop/333/all_diffsigCD8.csv",header =T)
m_deseqre_mo <- read.csv("C:/Users/a/Desktop/333/all_diffsigMOCD14.csv",header =T)
m_deseqre_nk <- read.csv("C:/Users/a/Desktop/333/all_diffsigNK.csv",header = T)

colnames(m_deseqre_cd4)[1] <- "hgnc_symbol"
colnames(m_deseqre_cd8)[1] <- "hgnc_symbol"
colnames(m_deseqre_mo)[1] <- "hgnc_symbol"
colnames(m_deseqre_nk)[1] <- "hgnc_symbol"

m_deseqre_cd4 <- merge(m_deseqre_cd4,t6g,by = "hgnc_symbol")
m_deseqre_cd8 <- merge(m_deseqre_cd8,t6g,by = "hgnc_symbol")
m_deseqre_mo <- merge(m_deseqre_mo,t6g,by = "hgnc_symbol")
m_deseqre_nk <- merge(m_deseqre_nk,t6g,by = "hgnc_symbol")


m_deseqsg_cd4 <- m_deseqre_cd4[abs(m_deseqre_cd4$padj)<0.001&abs(m_deseqre_cd4$log2FoldChange)>1,]
m_deseqsg_cd8 <- m_deseqre_cd8[abs(m_deseqre_cd8$padj)<0.001&abs(m_deseqre_cd8$log2FoldChange)>1,]
m_deseqsg_mo <- m_deseqre_mo[abs(m_deseqre_mo$padj)<0.001&abs(m_deseqre_mo$log2FoldChange)>1,]
m_deseqsg_nk <- m_deseqre_nk[abs(m_deseqre_nk$padj)<0.001&abs(m_deseqre_nk$log2FoldChange)>1,]

#使用正则表达式和grepl()函数检查列中的每个元素是否为数字或"X"和"Y"字符：
numeric_mask <- grepl("^[0-9XY]+$", m_deseqsg_cd4$chromosome_name)
m_deseqsg_cd4 <- subset(m_deseqsg_cd4, numeric_mask)
m_deseqsg_cd4$chromosome_name <- paste0("chr", m_deseqsg_cd4$chromosome_name)

numeric_mask <- grepl("^[0-9XY]+$", m_deseqsg_cd8$chromosome_name)
m_deseqsg_cd8 <- subset(m_deseqsg_cd8, numeric_mask)
m_deseqsg_cd8$chromosome_name <- paste0("chr", m_deseqsg_cd8$chromosome_name)

numeric_mask <- grepl("^[0-9XY]+$", m_deseqsg_mo$chromosome_name)
m_deseqsg_mo <- subset(m_deseqsg_mo, numeric_mask)
m_deseqsg_mo$chromosome_name <- paste0("chr", m_deseqsg_mo$chromosome_name)

numeric_mask <- grepl("^[0-9XY]+$", m_deseqsg_nk$chromosome_name)
m_deseqsg_nk <- subset(m_deseqsg_nk, numeric_mask)
m_deseqsg_nk$chromosome_name <- paste0("chr", m_deseqsg_nk$chromosome_name)



######将四种细胞的矩阵用subset函数分离出来，使用之前结果里的for循环
######5.18使用之前的grp做法


ascount <- read.csv("D:/data_anlysis/asncRNA/lnc_mrna_coexpression/ncrna_as_genename_count.csv",header  = T)
rownames(ascount) <- ascount$X
ascount$X <- NULL

#制造一个大致包含有全部lnc信息的charlist方便之后删除
lncallgene <- rownames(ascount)

m_deseq2sg_cd4<- m_deseqsg_cd4[!(m_deseqsg_cd4$hgnc_symbol %in%lncallgene ),]
m_deseq2sg_cd8<- m_deseqsg_cd8[!(m_deseqsg_cd8$hgnc_symbol %in%lncallgene ),]
m_deseq2sg_mo<- m_deseqsg_mo[!(m_deseqsg_mo$hgnc_symbol %in%lncallgene ),]
m_deseq2sg_nk<- m_deseqsg_nk[!(m_deseqsg_nk$hgnc_symbol %in%lncallgene ),]


#####对于count矩阵正式取子集
###########################################################################???下面的不要用
#先只用AS的样本

m_CD4countall <- cbind(m_CD4countas,m_CD4countc)
m_CD8countall <- cbind(m_CD8countas,m_CD8countc)
m_mocountall <- cbind(m_mocountas,m_mocountc)
m_nkcountall <- cbind(m_nkcountas,m_nkcountc)

lnc_CD4countall <- cbind(lnc_CD4countas,lnc_CD4countc)
lnc_CD8countall <- cbind(lnc_CD8countas,lnc_CD8countc)
lnc_mocountall <- cbind(lnc_mocountas,lnc_mocountc)
lnc_nkcountall <- cbind(lnc_nkcountas,lnc_nkcountc)


m_cd4_mat <- m_CD4countall[rownames(m_CD4countall) %in% m_deseq2sg_cd4, ]
m_cd8_mat <- m_CD8countall[rownames(m_CD8countall) %in% m_deseq2sg_cd8, ]
m_mo_mat <- m_mocountall[rownames(m_mocountall) %in% m_deseq2sg_mo, ]
m_nk_mat <- m_nkcountall[rownames(m_nkcountall) %in% m_deseq2sg_nk, ]

lnc_cd4_mat <- lnc_CD4countall[rownames(lnc_CD4countall) %in% lnc_deseqsg_cd4, ]
lnc_cd8_mat <- lnc_CD8countall[rownames(lnc_CD8countall) %in% lnc_deseqsg_cd8, ]
lnc_mo_mat <- lnc_mocountall[rownames(lnc_mocountall) %in% lnc_deseqsg_mo, ]
lnc_nk_mat <- lnc_nkcountall[rownames(lnc_nkcountall) %in% lnc_deseqsg_nk, ]




##############################################


library(circlize)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")
#initailize with customize chromosome track
set.seed(123)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)


#annotation for gene location
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genelist<-read.table("select.txt")
all_miRNAs <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol',
                                 'chromosome_name','start_position','end_position'), 
                    filters ='hgnc_symbol', values =genelist, mart = ensembl)
#over
circos.clear()
#load(system.file(package = "circlize", "extdata", "all_mRNA.csv
#                 "))
bed<-m_deseq2sg_cd4
#circos.genomicRainfall(bed)
#circos.genomicRainfall(bed, pch = 16, cex = 0.4, col = "#FF000080")
bed_lnc<-lnc_deseqsg_cd4
#circos.genomicRainfall(bed_lnc, pch = 16, cex = 0.4, col = "#0000FF80")

circos.track(factors = bed$chromosome_name, x = bed$start_position, y = -log10(bed$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "firebrick3")
             })
circos.track(factors = bed_lnc$chromosome_name, x = bed_lnc$start_position, y = -log10(bed_lnc$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "dodgerblue3")
             })
##circos.genomicPoints(region = ,value = ,)

#circos.track(factors = bed$chromosome_name, x = (bed$start_position+bed$end_position)/2, y = bed$log2FoldChange,
#             panel.fun = function(x, y) {
#               circos.points(x = x, y = y, pch = 16, cex = 0.5, col = "yellow")
#             })

# 检查 bed_lnc$log2FoldChange 是否有缺失值
any_na_log2FoldChange <- any(is.na(bed_lnc$padj))
print(any_na_log2FoldChange)




#######Cd8

library(circlize)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")

bed<-m_deseq2sg_cd8
#circos.genomicRainfall(bed)
#circos.genomicRainfall(bed, pch = 16, cex = 0.4, col = "#FF000080")
bed_lnc<-lnc_deseqsg_cd8
#circos.genomicRainfall(bed_lnc, pch = 16, cex = 0.4, col = "#0000FF80")

circos.track(factors = bed$chromosome_name, x = bed$start_position, y = -log10(bed$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "firebrick3")
             })
circos.track(factors = bed_lnc$chromosome_name, x = bed_lnc$start_position, y = -log10(bed_lnc$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "dodgerblue3")
             })
##circos.genomicPoints(region = ,value = ,)

circos.track(factors = bed$chromosome_name, x = (bed$start_position+bed$end_position)/2, y = log10(bed$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.5, col = "yellow")
             })




library(circlize)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")

bed<-m_deseq2sg_mo
#circos.genomicRainfall(bed)
#circos.genomicRainfall(bed, pch = 16, cex = 0.4, col = "#FF000080")
bed_lnc<-lnc_deseqsg_mo
#circos.genomicRainfall(bed_lnc, pch = 16, cex = 0.4, col = "#0000FF80")

circos.track(factors = bed$chromosome_name, x = bed$start_position, y = -log10(bed$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "firebrick3")
             })
circos.track(factors = bed_lnc$chromosome_name, x = bed_lnc$start_position, y = -log10(bed_lnc$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "dodgerblue3")
             })
##circos.genomicPoints(region = ,value = ,)

circos.track(factors = bed$chromosome_name, x = (bed$start_position+bed$end_position)/2, y = log10(bed$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.5, col = "yellow")
             })


library(circlize)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg19")

bed<-m_deseq2sg_nk
#circos.genomicRainfall(bed)
#circos.genomicRainfall(bed, pch = 16, cex = 0.4, col = "#FF000080")
bed_lnc<-lnc_deseqsg_nk
#circos.genomicRainfall(bed_lnc, pch = 16, cex = 0.4, col = "#0000FF80")

circos.track(factors = bed$chromosome_name, x = bed$start_position, y = -log10(bed$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "firebrick3")
             })
circos.track(factors = bed_lnc$chromosome_name, x = bed_lnc$start_position, y = -log10(bed_lnc$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.7, col = "dodgerblue3")
             })
##circos.genomicPoints(region = ,value = ,)

circos.track(factors = bed$chromosome_name, x = (bed$start_position+bed$end_position)/2, y = log10(bed$pvalue),
             panel.fun = function(x, y) {
               circos.points(x = x, y = y, pch = 16, cex = 0.5, col = "yellow")
             })


