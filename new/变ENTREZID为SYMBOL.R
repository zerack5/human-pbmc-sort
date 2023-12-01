library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(plyr)
library(dplyr)
library(parallel)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ReactomePA)

entrez_to_symbol = function(genes, ref){
  temp = strsplit(as.character(genes), "/")
  temp = lapply(temp, function(x) ref$SYMBOL[match(x, ref$ENTREZID)])
  temp = lapply(temp, function(x) paste(x, collapse = "/"))
  temp = as.vector(t(as.data.frame(temp)))
  return(temp)
}

filename <- "D:/extrawork/Chondrocyte_2_GOpathway.xlsx"
filename <- "D:/extrawork/Chondrocyte_4_GOpathway.xlsx"
filename <- "D:/extrawork/Fibroblast_GOpathway.xlsx"


library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')


t2g <- biomaRt::getBM(attributes = c("entrezgene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, ENTREZID = entrezgene_id, SYMBOL = external_gene_name)
#colnames(t2g) = c("ENTREZID","SYMBOL")


sample_data <- readxl::read_xlsx(filename,sheet = "GO")

sample_data$core_enrichment <- entrez_to_symbol(sample_data$core_enrichment, ref = t2g)
write.xlsx(sample_data,filename)






sample_data <- readxl::read_xlsx(filename,sheet = )

sample_data$core_enrichment <- entrez_to_symbol(sample_data$core_enrichment, ref = t2g)
write.xlsx(sample_data,filename)


