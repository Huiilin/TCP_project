rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 200)

load("6_Multi_omics_enrich_ranking/forPlot_plasma.RData")

library(ggplot2)
library(tidyverse)
library(ggthemes)
library(export)
library(readxl)

### read in pathway class files
path = "6_Multi_omics_enrich_ranking/ReactomePathwaysRelation_downloaded.xlsx"
pathclass <- read_excel(path = path,sheet = 1,col_names = T)
colnames(pathclass) <- c("class","PathwayIdentifier")

### combine pathway class
names(plasma4)
names(plasma4)[which(names(plasma4) == "Pathway identifier")]  <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma4,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier0"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier1"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier2"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier3"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier4"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier5"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier6"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

names(plasma5)[which(names(plasma5) == "PathwayIdentifier")] <- "PathwayIdentifier7"
names(plasma5)[which(names(plasma5) == "class")] <- "PathwayIdentifier"
plasma5 <- right_join(pathclass,plasma5,by = "PathwayIdentifier")

a = grep(pattern = "PathwayIdentifier",colnames(plasma5))
plasma6 <- plasma5 %>% as.data.frame()

for (j in 1:nrow(plasma6)) {
  for (i in a) {
    if (is.na(plasma6$class[j]) == FALSE)
      next;
    if (is.na(plasma6[j,i]) == FALSE) 
      plasma6$class[j] <- plasma6[j,i]
  }
}

table(plasma6$class)

plasma6$pcName[plasma6$class == "R-HSA-1430728"] <- 'Metabolism'
plasma6$pcName[plasma6$class == "R-HSA-162582"] <- 'Signal Transduction'
plasma6$pcName[plasma6$class == "R-HSA-168256"] <- "Immune System"
plasma6$pcName[plasma6$class == "R-HSA-382551"] <- 'Transport of small molecules'
plasma6$pcName[plasma6$class == "R-HSA-392499"] <- 'Metabolism of proteins'
plasma6$pcName[plasma6$class == "R-HSA-5653656"] <- "Vesicle-mediated transport"
plasma6$pcName[plasma6$class == "R-HSA-8963743"] <- "Digestion and absorption"

table(plasma6$contrast)
plasma6 <- filter(plasma6,contrast != "2vs1")

### plot #######################################################################
plasma6 %>% dplyr::select(-(PathwayIdentifier:PathwayIdentifier1)) -> plasma7
plasma7 <- arrange(plasma7,pcName,PathwayIdentifier0)

which(plasma7$`Pathway name` == "Metabolism")
plasma8 <- plasma7[-4,]

plasma8 <- filter(plasma8,`Entities pValue` < 0.01)

a1 = table(plasma8$PathwayIdentifier0) %>% as.data.frame()

plasma8$IDcode <- rep(0,nrow(plasma8))
a1$index <- rep(0,nrow(a1))

for (i in 1:nrow(a1)) {
  for (j in 1:nrow(plasma8)) {
    if (plasma8$PathwayIdentifier0[j] == a1[i,1])
      a1[i,3] <- j
  }
}

a1 <- arrange(a1,index)
a1$index2 <- seq(1:nrow(a1))

for (j in 1:nrow(plasma8)) {
  for (i in 1:nrow(a1)) {
    if (plasma8$PathwayIdentifier0[j] == a1$Var1[i])
      plasma8$IDcode[j] <- a1$index2[i]
  }
}

plasma8[,c("Pathway name","IDcode")] %>% arrange(.,IDcode)
plasma8[,c("Pathway name","pcName","IDcode")] %>% arrange(.,pcName) -> IDcode

## plot
IDcode <- IDcode[!duplicated(IDcode$IDcode),]
plasma8$logP <- -log10(plasma8$`Entities pValue`)


## data rearrange -- median 
library(MASS)
xP = plasma8$median[plasma8$median > 0]
yP = 20000/3400 * xP + 3359/3400

plasma8$medianChanged[plasma8$median > 0]  <- 
  20000/3400 * plasma8$median[plasma8$median > 0] + 3359/3400

xN = plasma8$median[plasma8$median < 0]
yN = 80/219*xN -215/219

plasma8$medianChanged[plasma8$median < 0] <- 
  80/219*plasma8$median[plasma8$median < 0] -215/219

## data rearrange -- absweight
plasma8$weightChange = 0.781*plasma8$absWeightSums + 0.099

## pic with label 
pPLT1 <- ggplot(plasma8,aes(x = contrast,y = IDcode)) +
  geom_point(aes(size = weightChange,colour = logP,fill = medianChanged),
             shape =21,stroke = 1) +
  scale_y_continuous(breaks = 1:max(plasma8$IDcode),
                     labels = IDcode$`Pathway name`,
                     sec.axis = sec_axis(~.,breaks = 1:max(plasma8$IDcode),
                                         labels = IDcode$pcName),na.value = "") +
  scale_fill_gradient2(low = '#3288bd',high = '#d53e4f',name = "Pathway Direction")+
  scale_color_gradient2(low = "#e0e0e0",high = "black",name = "-log(P)") +
  scale_size_area(max_size = 10,name = "Pathway Weight") +
  ggtitle("Pathway Ranking of Plasma") +
  theme_few(base_size = 18) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(5,"pt"),
        plot.margin = margin(c(0,0,0,0)))
pPLT1

## pic without label
pPLT2 <- ggplot(plasma8,aes(x = contrast,y = IDcode)) +
  geom_point(aes(size = weightChange,colour = logP,fill = medianChanged),
             shape =21,stroke = 1) +
  scale_y_continuous(breaks = 1:max(plasma8$IDcode),
                     sec.axis = sec_axis(~.,breaks = 1:max(plasma8$IDcode)),
                     na.value = "") +
  scale_fill_gradient2(low = '#3288bd',high = '#d53e4f',name = "Pathway Direction",
                       midpoint = 0)+
  scale_color_gradient2(low = "#e0e0e0",high = "black",name = "-log(P)") +
  scale_size_area(max_size = 10,name = "Pathway Weight") +
  ggtitle("Pathway Ranking of Plasma") +
  theme_few(base_size = 18) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(5,"pt"),
        plot.margin = margin(c(0,0,0,0)))
pPLT2

## save
save.image(file = "6_Multi_omics_enrich_ranking/figures_v2_plasma.RData")
graph2ppt(pPLT2,file = "0_figures/figures_v3-0.pptx",append = TRUE,width = 6.5,height= 20,
          bg = "transparent")
graph2ppt(pPLT1,file = "0_figures/figures_v3-0.pptx",append = TRUE,width = 15,height= 15,
          bg = "transparent")
