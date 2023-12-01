
setwd("C:/Users/a/Desktop/表达矩阵")

library(DESeq2)
library(dplyr)
#导入counts数据矩阵，以行为基因，列为样本

filename1 <- "D:/bulk_control_AS/svm/cd4_count.csv"
data2as <- read.csv('expas.csv', sep = ',',header = T,row.names = 1)
sdata2as <- data2as[c(grep("CD4",colnames(data2as)))]

data2c <- read.csv('expcontrol.csv', sep = ',',header = T,row.names = 1)
sdata2c <- data2c[c(grep("CD4",colnames(data2c)))]

scb2data <- cbind(sdata2as, sdata2c)

write.csv(scb2data,filename1)

######################################################cd8
filename1 <- "D:/bulk_control_AS/svm/cd8_count.csv"

data2as <- read.csv('expas.csv', sep = ',',header = T,row.names = 1)
sdata2as <- data2as[c(grep("CD8",colnames(data2as)))]

data2c <- read.csv('expcontrol.csv', sep = ',',header = T,row.names = 1)
sdata2c <- data2c[c(grep("CD8",colnames(data2c)))]


scb2data <- cbind(sdata2as, sdata2c)

write.csv(scb2data,filename1)

################################################################mo

filename1 <- "D:/bulk_control_AS/svm/mo_count.csv"

data2as <- read.csv('expas.csv', sep = ',',header = T,row.names = 1)
sdata2as <- data2as[c(grep("_M",colnames(data2as)))]

data2c <- read.csv('expcontrol.csv', sep = ',',header = T,row.names = 1)
sdata2c <- data2c[c(grep("_M",colnames(data2c)))]

scb2data <- cbind(sdata2as, sdata2c)

write.csv(scb2data,filename1)


#######################################nk

filename1 <- "D:/bulk_control_AS/svm/nk_count.csv"
data2as <- read.csv('expas.csv', sep = ',',header = T,row.names = 1)
sdata2as <- data2as[c(grep("_N",colnames(data2as)))]

data2c <- read.csv('expcontrol.csv', sep = ',',header = T,row.names = 1)
sdata2c <- data2c[c(grep("_N",colnames(data2c)))]

scb2data <- cbind(sdata2as, sdata2c)

write.csv(scb2data,filename1)





