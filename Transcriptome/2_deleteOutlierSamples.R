rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

library(tidyverse)
library(data.table)
library(readxl)

load("3_RNA_V2/Trans_medinfo.RData")

## read in FPKM data
file = "3_RNA_V2/gene_sample_FPKM.xlsx"
FPKM <- readxl::read_excel(path = file, sheet = 1,col_names = TRUE) %>% as.data.frame
rownames(FPKM) <- unlist(FPKM$geneid)

save(FPKM, file = "3_RNA_V2/gene_FPKM.RData")

FPKM <- FPKM[,-1]
FPKM <- FPKM[,rownames(medInfoTrans)]
names(FPKM)

## remove outlier samples 
which(names(FPKM)=="A2") # 55
which(names(FPKM)=="D3") # 41
which(names(FPKM)=="F1") # 56

FPKM54 <- FPKM[,c(-55,-41,-56)]

which(rownames(medInfoTrans) == "A2") # 55
which(rownames(medInfoTrans) == "D3") # 41
which(rownames(medInfoTrans) == "F1") # 56

medInfoTrans <- medInfoTrans[c(-55,-41,-56),]

all(rownames(medInfoTrans) == colnames(FPKM54))

save(FPKM54,medInfoTrans,file = "3_RNA_V2/gene_FPKM_54.RData")
