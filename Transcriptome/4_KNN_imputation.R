rm(list = ls())
options(scipen = 200)

load("gene_FPKM_54_outliers2na.RData")

library(impute)
library(tidyverse)

## impute.KNN requires input data samples in columns, genes in rows
## our input data: sample in rows 
## so need to t(data) first
cols = str_subset(colnames(RNA_outliers2na),pattern = "ENSG")

imputeDat <- t(RNA_outliers2na[,cols]) %>% as.data.frame

imputeDat1 <- apply(imputeDat,2,as.numeric)  %>% as.data.frame
rownames(imputeDat1) = rownames(imputeDat)

## impute
aa = impute.knn(as.matrix(imputeDat1), maxp = 50000)
imputeDat2 <- aa$data

table(is.na(imputeDat2))

save(imputeDat2,file = "gene_FPKM_54_KNNimputed.RData")
