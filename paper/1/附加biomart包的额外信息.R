

path0 <- "D:/data_anlysis/asncRNA"
setwd(path0)


file1 <- "1phase_il17acount.csv"
file2 <- "2phase_il17acount.csv"
file3 <- "3phase_il17acount.csv"

de_cd4 <- "deseqresult/deseq2_CD4.csv"
de_cd4_sig <- "deseqresult/deseq2_CD4_sig.csv"
de_cd4_sig_up <- "deseqresult/deseq2_CD4_sig_up.csv"
de_cd4_sig_down <- "deseqresult/deseq2_CD4_sig_down.csv"

de_cd8 <- "deseqresult/deseq2_CD8.csv"
de_cd8_sig <- "deseqresult/deseq2_CD8_sig.csv"
de_cd8_sig_up <- "deseqresult/deseq2_CD8_sig_up.csv"
de_cd8_sig_down <- "deseqresult/deseq2_CD8_sig_down.csv"

de_mo <- "deseqresult/deseq2_mo.csv"
de_mo_sig <- "deseqresult/deseq2_mo_sig.csv"
de_mo_sig_up <- "deseqresult/deseq2_mo_sig_up.csv"
de_mo_sig_down <- "deseqresult/deseq2_mo_sig_down.csv"

de_nk <- "deseqresult/deseq2_nk.csv"
de_nk_sig <- "deseqresult/deseq2_nk_sig.csv"
de_nk_sig_up <- "deseqresult/deseq2_nk_sig_up.csv"
de_nk_sig_down <- "deseqresult/deseq2_nk_sig_down.csv"


####cd4
diff<-read.csv(file=de_cd4)
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
t5g <- readRDS("D:/Rdata/biomart/t5g.rdata")
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
write.csv(dataf,de_cd4,quote = F,row.names = F,col.names = T)

######cd8
diff<-read.csv(file=de_cd8)
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
t5g <- readRDS("D:/Rdata/biomart/t5g.rdata")
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
write.csv(dataf,de_cd8,quote = F,row.names = F,col.names = T)

######mo
diff<-read.csv(file=de_mo)
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
t5g <- readRDS("D:/Rdata/biomart/t5g.rdata")
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
write.csv(dataf,de_mo,quote = F,row.names = F,col.names = T)

######nk
diff<-read.csv(file=de_nk)
library(dplyr)
diff <- arrange(diff,padj)
colnames(diff)[1] <- "ensembl_gene_id"
t5g <- readRDS("D:/Rdata/biomart/t5g.rdata")
tx2g <- t5g[,c("ensembl_gene_id","external_gene_name","entrezgene_id")]
tx2g <- unique(tx2g)
dataf <- merge(diff,tx2g,by="ensembl_gene_id")
dataf <- arrange(dataf,desc(log2FoldChange))
write.csv(dataf,de_nk,quote = F,row.names = F,col.names = T)





