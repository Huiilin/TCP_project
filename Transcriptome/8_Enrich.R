rm(list = ls())
options(scipen = 200)

load("5_MOFA2_V5/RData_FCfilter/RNA_All.RData")

library(ReactomePA)
library(org.Hs.eg.db)
library(tidyverse)

## ID convert
ENSG <- keys(org.Hs.eg.db, keytype = "ENSEMBL")

list <- AnnotationDbi::select(org.Hs.eg.db,keys = ENSG,
                              columns = c("ENTREZID","SYMBOL"),
                              keytype = "ENSEMBL")

## Enrich function
enrich <- function(ID) {
  IDList <- list[match(ID,list[,"ENSEMBL"]),]
  IDList <- na.omit(IDList)
  ENTREZ <- IDList$ENTREZID
  
  Enrich <- enrichPathway(ENTREZ,
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1)
  pathway <- Enrich@result
  pathway <- dplyr::filter(pathway,pvalue < 0.05)
  
  return(pathway)
}

## RNA
data <- list(diff10,diff20,diff21)
pathway <- list()

for (i in 1:3) {
  a <- as.data.frame(data[i])
  ID <- a$featureID
  b <- enrich(ID)
  pathway[i] <- list(b)
}

## writeFiles
rm(a,b,i,ID)

tmp <- paste(getwd(),"5_MOFA2_V5/Enrich/",sep = "/")
tmp
pathway10 <- pathway[[1]] %>% mutate(., contrast = rep("0 vs 1",nrow(.)))
pathway20 <- pathway[[2]] %>% mutate(., contrast = rep("0 vs 2",nrow(.)))
pathway21 <- pathway[[3]] %>% mutate(., contrast = rep("1 vs 2",nrow(.)))

pathAll <- rbind(pathway10,pathway20) %>% rbind(.,pathway21)

## writeFiles
write.csv(pathAll, file = paste(tmp,"EnrichRes_RNA.csv", sep = ""), row.names = F)

## saveData
save.image(file = "5_MOFA2_V5/Enrich/Enrich_RNA.RData")
