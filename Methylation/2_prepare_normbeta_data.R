rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

library(tidyverse)

load("4_Methyl/ME_medinfo.RData")

## read in beta value
file = "4_Methyl/NormBeta_20191217.csv"
beta <- read_csv(file = file,col_names = TRUE)  %>% as.data.frame
names(beta)[1] <- "probeID"
rownames(beta) <- beta$probeID
beta <- beta[,-1]

## combine samples
a <- intersect(names(beta),medinfoME$sampleID)
beta56 <- select(beta, all_of(a))
medinfoME56 <- medinfoME[which(medinfoME$sampleID %in% a),]
beta56 <- beta56[,medinfoME56$sampleID]

all(medinfoME56$sampleID %in% names(beta56))
all(medinfoME56$sampleID == names(beta56))

save(beta56,medinfoME56,file = "4_Methyl/NormBeta56.RData")
