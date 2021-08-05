rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

library(tidyverse)

load("gene_FPKM_54.RData")
load("IQRNumFunc_RNA.RData")
load("IQROutlierFunc_RNA.RData")
load("featureSelect_Func.RData")
load("outliers2NAFunc.RData")

## remove low FPKM
keep <- rowSums(FPKM54 == 0) == 54
table(keep)

FPKM54 <- FPKM54[!keep,]

## logNorm
logNorm <- log2(FPKM54 + 1)

save(logNorm, file = "3_RNA_V2/gene_FPKM_54_logNorm.RData")

dataRNA <- t(logNorm) %>% as.data.frame
dataRNA$sampleID <- rownames(dataRNA)

## merge FPKM with clinical data
G <- medInfoTrans[,c(14,15)]
RNAGroup <- merge(G,dataRNA,by = "sampleID")

## seperate groups
RNAGroup0 <- RNAGroup[which(RNAGroup$TCP == 0),]
RNAGroup1 <- RNAGroup[which(RNAGroup$TCP == 1),]
RNAGroup2 <- RNAGroup[which(RNAGroup$TCP == 2),]

## calclulate IQR for each variables
RNAIQR0 <- IQRNum(RNAGroup0)
RNAIQR1 <- IQRNum(RNAGroup1)
RNAIQR2 <- IQRNum(RNAGroup2)

## set outliers based on IQR value
RNAoutliers0 <- outliers(RNAIQR0)
RNAoutliers1 <- outliers(RNAIQR1)
RNAoutliers2 <- outliers(RNAIQR2)

RNAoutliers0 <- RNAoutliers0[1:39,]
RNAoutliers1 <- RNAoutliers1[1:9,]
RNAoutliers2 <- RNAoutliers2[1:6,]

## select outliers features
all <- rbind(RNAoutliers0,RNAoutliers1) %>% rbind(.,RNAoutliers2)
RNA_outliers2na <- featureSelect(all)

keep <- colnames(RNA_outliers2na)[which(RNA_outliers2na[55,] == "keep")]
keep <- append(keep,c("sampleID" ,"TCP" ),0)

RNA_outliers2na <- RNA_outliers2na[,keep]

## outliers to NA
RNA_outliers2na <- outliers2NA(RNA_outliers2na)
rownames(RNA_outliers2na) <- RNA_outliers2na$sampleID

save(RNA_outliers2na, file = "gene_FPKM_54_outliers2na.RData")
