library(tidyverse)

## input data: resdata is a list of res10,res20,res21
##             deData is limma/champ analysed DE data

FCfilter <- function(resdata,diffMetabos10,diffMetabos20,diffMetabos21,deData,IDAll) {

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
  
  names(FCAll)[2] <- "logFC01"
  names(FCAll)[3] <- "logFC02"
  names(FCAll)[4] <- "logFC12"
  
  FCAll_1 <- dplyr::filter(FCAll,logFC01 > 0, logFC12 > 0)
  FCAll_2 <- dplyr::filter(FCAll,logFC01 < 0, logFC12 < 0)
  
  FCAll_3 <- rbind(FCAll_1,FCAll_2)
  
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

  keep <- FCAll_3$featureID
  
  deMofa <- deData[keep,]
  
  delist <- list("diff10" = diff10, "diff20" = diff20, "diff21" = diff21,
                 "diffAll" = diff,"deData" = deMofa,"filter" = FCAll_3)
  
  return(delist)
}
