rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 200)

library(tidyverse)
library(matrixStats)

## transform data 
dataME <- t(beta56) %>% as.data.frame
dataME$sampleID <- rownames(dataME)

## merge
G <- medinfoME56[,c(14,15)]
MEGroup <- merge(G,dataME,by = "sampleID")

MEGroup0 <- MEGroup[which(MEGroup$TCP ==0),]
MEGroup1 <- MEGroup[which(MEGroup$TCP ==1),]
MEGroup2 <- MEGroup[which(MEGroup$TCP ==2),]

rownames(MEGroup0) <- MEGroup0$sampleID
rownames(MEGroup1) <- MEGroup1$sampleID
rownames(MEGroup2) <- MEGroup2$sampleID

## Group0 CV
MEGroup0 <- MEGroup0[,c(-1,-2)]
MEGroup00 <- t(MEGroup0) %>% as.data.frame
head(MEGroup0)[1:2,1:6]
head(MEGroup00)[1:2,1:6]
MEGroup00$probeID <- rownames(MEGroup00)

cols <- names(MEGroup00)[1:39]

QCdata0 <- MEGroup00 %>% filter(rowSums(.[cols]) != 0) %>% 
  mutate(mean = rowMeans(.[cols]), sd = rowSds(as.matrix(.[cols]))) %>%
  mutate(CV = sd/mean)

qc_sig01 <- QCdata0 %>% filter(.["CV"] < 0.3)

## Group1 CV
MEGroup1 <- MEGroup1[,c(-1,-2)]
MEGroup10 <- t(MEGroup1) %>% as.data.frame
head(MEGroup10)
MEGroup10$probeID <- rownames(MEGroup10)

cols <- names(MEGroup10)[1:10]

QCdata1 <- MEGroup10 %>% filter(rowSums(.[cols]) != 0) %>% 
  mutate(mean = rowMeans(.[cols]), sd = rowSds(as.matrix(.[cols]))) %>%
  mutate(CV = sd/mean)
 
qc_sig11 <- QCdata1 %>% filter(.["CV"] < 0.3)

## Group2 CV
MEGroup2 <- MEGroup2[,c(-1,-2)]
MEGroup20 <- t(MEGroup2) %>% as.data.frame
head(MEGroup20)[1:2,1:6]

MEGroup20$probeID <- rownames(MEGroup20)

cols <- names(MEGroup20)[1:7]

QCdata2 <- MEGroup20 %>% filter(rowSums(.[cols]) != 0) %>% 
  mutate(mean = rowMeans(.[cols]), sd = rowSds(as.matrix(.[cols]))) %>%
  mutate(CV = sd/mean)
  
qc_sig21 <- QCdata2 %>% filter(.["CV"] < 0.3)

## CV < 0.3
ID01 <- qc_sig01$probeID
ID11 <- qc_sig11$probeID
ID21 <- qc_sig21$probeID

IDALL2 <- intersect(ID01,ID11) %>% intersect(.,ID21)
beta56CV2 <- beta56[IDALL2,]

## save
save(beta56CV2,file = "Methyl_QC_CV_03_intersect.RData")



