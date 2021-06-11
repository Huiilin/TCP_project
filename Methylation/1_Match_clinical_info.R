rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

library(tidyverse)

load("1_clinical_info/medInfo_covariatesFactorised.RData")
load("sampleInfoAll.RData")

# methylation 
sampleInfoME <- sampleInfo[,c(1,5,10)]
str(sampleInfoME)
names(sampleInfoME) <- c("medicalID","sampleID","TCP")

medinfoME <- merge(medInfo,sampleInfoME, by = "medicalID")

a = which(names(medinfoME) == "PatientID")
b = which(names(medinfoME) == "Name")

medinfoME <- medinfoME[,c(-a,-b)]

medinfoME <- medinfoME[complete.cases(medinfoME),] %>% .[!duplicated(.),]

medinfoME$condition[medinfoME$TCP == 0] <- "NoTCP"
medinfoME$condition[medinfoME$TCP == 1] <- "TCP1"
medinfoME$condition[medinfoME$TCP == 2] <- "TCP2"

medinfoME$condition <- as.factor(medinfoME$condition)
rownames(medinfoME) <- unlist(medinfoME$sampleID)

save(medinfoME,sampleInfoME, file = "4_Methyl/ME_medinfo.RData")
