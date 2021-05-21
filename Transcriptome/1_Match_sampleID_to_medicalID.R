rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

library(tidyverse)

load("1_clinical_info/medInfo_covariatesFactorised.RData")
load("sampleInfoAll.RData")

## transcriptome ###
sampleInfoTrans <- sampleInfo[,c(1,6,10)]
str(sampleInfoTrans)
names(sampleInfoTrans) <- c("medicalID","sampleID","TCP")

medInfoTrans <- merge(medInfo,sampleInfoTrans,by = "medicalID")

a = which(names(medInfoTrans) == "PatientID")
b = which(names(medInfoTrans) == "Name")

medInfoTrans <- medInfoTrans[,c(-a,-b)]

medInfoTrans <- medInfoTrans[complete.cases(medInfoTrans),] %>% .[!duplicated(.),]
  
medInfoTrans$condition[medInfoTrans$TCP == 0] <- "NoTCP"
medInfoTrans$condition[medInfoTrans$TCP == 1] <- "TCP1"
medInfoTrans$condition[medInfoTrans$TCP == 2] <- "TCP2"
medInfoTrans$condition <- as.factor(medInfoTrans$condition)
rownames(medInfoTrans) <- unlist(medInfoTrans$sampleID)

save(medInfoTrans,sampleInfoTrans,file = "3_RNA_V2/Trans_medinfo.RData")  
