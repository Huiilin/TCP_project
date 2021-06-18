rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

library(tidyverse)
library(ChAMP)

## load in data
load("4_Methyl/NormBeta56.RData")
load("4_Methyl/Methyl_QC_CV_03_intersect.RData")

names(medinfoME56)[which(names(medinfoME56) == "sampleID")] <- "Sample_Name"
names(medinfoME56)[which(names(medinfoME56) == "condition")] <- "Sample_Group"

## create methylation list
beta56CV2 <- as.matrix(beta56CV2)
medinfoME56$Sample_Group <- as.character(medinfoME56$Sample_Group)
names(medinfoME56)
medinfoME56Ch <- medinfoME56[,c(1,2,4,9:14,16)]

## CHAMP
mydata <- list(beta = beta56CV2, pd = medinfoME56Ch)

myCombat <- champ.runCombat(mydata$beta,
                            pd = mydata$pd,
                            variablename = "Sample_Group",
                            batchname = c("HistologyType","Sex","ageGroup",
                                          "HeightGroup","WeightGroup",
                                          "Staging"))
                                         
## differential analysis
myDMP <- champ.DMP(myCombat,pheno = mydata$pd$Sample_Group,arraytype = "EPIC",adjPVal = 1)

res10 <- myDMP$NoTCP_to_TCP1
diff10 <- res10[which(abs(res10$deltaBeta) > 0.1 & res10$P.Value < 0.05),]

res20 <- myDMP$NoTCP_to_TCP2
diff20 <- res20[which(abs(res20$deltaBeta) > 0.1 & res20$P.Value < 0.05),]

res21 <- myDMP$TCP1_to_TCP2
diff21 <- res21[which(abs(res21$deltaBeta) > 0.1 & res21$P.Value < 0.05),]

## linear FC
source("5_MOFA2_V5/0_newcutoff_linearFC_fun.R")                  

select(diff10,-c("NoTCP_AVG",'TCP1_AVG')) -> diff10
select(diff20,-c("NoTCP_AVG",'TCP2_AVG')) -> diff20
select(diff21,-c("TCP1_AVG",'TCP2_AVG')) -> diff21

resdata <- list("0vs1" = res10,"0vs2" = res20, "1vs2" = res21)
Mlist <- FCfilter(resdata,diff10,diff20,diff21,beta_DE,IDAll)

FCAll <- Mlist$filter

diff <- Mlist$diffAll
diff <- diff[,c(20,1:7,19,8:18)]

Methyl_MOFA <- Mlist$deData

## save
save(Methyl_MOFA,medinfoME56,file = "Methyl_FCfilter_DE.RData")
save(Mlist,file = "Methyl_All.RData")

## PCA
source("6_PCA_fun.R")

names(medinfoME56)[which(names(medinfoME56) == "Sample_Name")] <- "sampleID"

p <- metaboPCA(Methyl_MOFA,medinfoME56)
p

ggsave(filename = "5_MOFA2_V5/PCA_Image/Methyl_PCA.png",width = 5,height = 5)

