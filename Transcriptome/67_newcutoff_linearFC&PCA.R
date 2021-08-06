rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 200)

library(tidyverse)
##################################################################
### select FC increase or decrease features among 3 groups ######
##################################################################
load("RNA_limma_1219.RData")

resdata <- list("0vs1" = res10,"0vs2" = res20, "1vs2" = res21)
FCdata <- list()

for (i in 1:3) {
  a = as.data.frame(resdata[i])
  a <- mutate(a, featureID = rownames(a))
  resdata[[i]] <-a
  FCdata[i] <- list(dplyr::filter(a, featureID %in% IDAll))
}

names(FCdata) <- c("0vs1","0vs2","1vs2")

FC01 <- FCdata[[1]]
FC02 <- FCdata[[2]]
FC12 <- FCdata[[3]]

FCAll <- merge(FC01,FC02,by = "featureID") %>% 
  merge(.,FC12,by ="featureID" ) %>% 
  dplyr::select(.,c("featureID","X0vs1.logFC","X0vs2.logFC","X1vs2.logFC"))

names(FCAll)
names(FCAll)[2] <- "logFC01"
names(FCAll)[3] <- "logFC02"
names(FCAll)[4] <- "logFC12"

FCAll_1 <- dplyr::filter(FCAll,logFC01 > 0, logFC12 > 0)
FCAll_2 <- dplyr::filter(FCAll,logFC01 < 0, logFC12 < 0)

FCAll_3 <- rbind(FCAll_1,FCAll_2)

## filter group results
diff10 <- dplyr::filter(diffMetabos10,rownames(diffMetabos10) %in% FCAll_3$featureID)
diff20 <- dplyr::filter(diffMetabos20,rownames(diffMetabos20) %in% FCAll_3$featureID)
diff21 <- dplyr::filter(diffMetabos21,rownames(diffMetabos21) %in% FCAll_3$featureID)

diff10$contrast <- rep("0 vs 1",nrow(diff10))
diff20$contrast <- rep("0 vs 2",nrow(diff20))
diff21$contrast <- rep("1 vs 2",nrow(diff21))

diff10$featureID <- rownames(diff10)
diff20$featureID <- rownames(diff20)
diff21$featureID <- rownames(diff21)

diff <- rbind(diff10,diff20) %>% rbind(.,diff21)
diff <- diff[,c(8,1:6,7)]
keep <- FCAll_3$featureID
RNA_MOFA <- RNA_DE[keep,]

## saveFiles
save(RNA_MOFA,RNA,file = "RNA_FCfilter_DE.RData")

## PCA 
source("6_PCA_fun.R")
load("Trans_medinfo.RData")

RNA <- dplyr::filter(medInfoTrans,sampleID %in% colnames(RNA_MOFA))
pRNA <- metaboPCA(RNA_MOFA,RNA)
