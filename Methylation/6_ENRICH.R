rm(list = ls())
options(scipen = 200)

load("Methyl_All.RData")

library(ReactomePA)
library(org.Hs.eg.db)
library(tidyverse)

ENSG <- keys(org.Hs.eg.db, keytype = "ENSEMBL")

list <- AnnotationDbi::select(org.Hs.eg.db,keys = ENSG,
                              columns = c("ENTREZID","SYMBOL"),
                              keytype = "ENSEMBL")
                              
## enrich functin for methyl
enrich <- function(ID) {
  IDList <- list[match(ID,list[,"SYMBOL"]),]
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

## get genes
diffM10 <- Mlist$diff10
diffM20 <- Mlist$diff20
diffM21 <- Mlist$diff21

data <- list(diffM10,diffM20,diffM21)
pathway <- list()

for (i in 1:3) {
  a <- as.data.frame(data[i])
  gene <- a$gene
  pathway[i] <- list(enrich(gene))
}

## writeFiles

tmp <- paste(getwd(),"path",sep = "/")

pathway10 <- pathway[[1]] %>% mutate(., contrast = rep("0 vs 1",nrow(.)))
pathway20 <- pathway[[2]] %>% mutate(., contrast = rep("0 vs 2",nrow(.)))
pathway21 <- pathway[[3]] %>% mutate(., contrast = rep("1 vs 2",nrow(.)))

pathAll <- rbind(pathway10,pathway20) %>% rbind(.,pathway21)

## writeFiles
write.csv(pathAll, file = paste(tmp,"EnrichRes_Methyl.csv", sep = ""), 
          row.names = F)
          
## saveData
save.image(file = "Enrich_Methyl.RData")
