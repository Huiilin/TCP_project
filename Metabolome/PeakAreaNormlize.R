rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 200)

library(readxl)
library(tidyverse)

path2 = "inputFile_V2/LMM-定性KEGG&HMDB物质查找-含甲基化硫化氧化标注_20201222.xlsx"
sheet = "combined_selected"
cmpd = read_excel(path = path2, sheet = sheet,col_names = T)

####################
#####         ######
##### Tissue  ######
#####         ######
####################

## pos
path = "inputFile_V2/CRT_pos_IS_adj.xlsx"
TissueP <- read_excel(path = path,sheet = 1,col_names = T)
TissueP <- select(TissueP,-(Compounds:Name))

## neg
path1 <- "inputFile_V2/CRT_neg_IS_adj.xlsx"
TissueN <- read_excel(path = path1,sheet = 1, col_names = T)
TissueN <- select(TissueN,-(Compounds:name))

## new medinfo tissue
load("tissue_medinfo.RData")
a = which(medInfoTissue$sampleID == "L12")

medInfoTissue[59,] <- medInfoTissue[a,]
medInfoTissue$sampleID[59] <- "L12_Pink"
medInfoTissue$sampleID[a] <- "L12_Yellow"

which(medInfoTissue$medicalID == "741044")
medInfoTissue$medicalID[18] <- "741044_1"

save(medInfoTissue,file = "tissue_medinfo_1.RData")

## select data 
TissueP1 <- select(TissueP,-(QC1:QC9)) %>% 
  select(.,UniqueID,contains(medInfoTissue$sampleID))

TissueN1 <- select(TissueN,-(QC1_20200716183204:QC9)) %>%
  select(.,UniqueID,contains(medInfoTissue$sampleID))

## new peak area normalization matrix 
## pos
rm(path,path1,path2)
n = dim(TissueP1)[2]

PA_TP <- matrix(0,nrow = dim(TissueP1)[1],ncol = (n-1),
                dimnames = list(TissueP1$UniqueID,names(TissueP1)[2:n])) %>% 
  as.data.frame

for (j in 1:(n-1)) {
  for (i in 1:dim(TissueP1)[1]) {
    PA_TP[i,j] <- TissueP1[i,(j+1)]/sum(TissueP1[,(j+1)])
  }
}

PA_TP$uniqueID <- rownames(PA_TP)

## neg
rm(i,j,n)
n = dim(TissueN1)[2]

PA_TN <- matrix(0,nrow = dim(TissueN1)[1], ncol = (n-1),
                dimnames = list(TissueN1$UniqueID,names(TissueN1)[2:n])) %>%
  as.data.frame

for (j in 1:(n-1)) {
  for (i in 1:dim(TissueN1)[1]) {
    PA_TN[i,j] <- TissueN1[i,(j+1)]/sum(TissueN1[,(j+1)])
  }
}

PA_TN$uniqueID <- rownames(PA_TN)

### remove NaN and Inf
a = NA

for (i in 1:(n-1)) {
  if (PA_TN[1,i] == "NaN") {
    a = append(a,TissueN1$UniqueID[which(TissueN1[,i] == "NaN")],0)
  }
}

TissueN1 <- filter(TissueN1,!(UniqueID %in% a))

## pick repeated chemicals
rm(i,j,n)
cmpdP <- filter(cmpd,Tissue == "CRT",adduct_type == "pos")
cmpdN <- filter(cmpd,Tissue == "CRT",adduct_type == "neg")

keepTP <- NA
keepTN <- NA

for (i in 1:dim(cmpdP)[1]) {
  for (j in 1:dim(cmpdN)[1]) {
    if (cmpdP$`English Name`[i] == cmpdN$`English Name`[j]) {
      keepTP <- append(keepTP,cmpdP$UniqueID[i],0)
      keepTN <- append(keepTN,cmpdN$UniqueID[j],0)
    }
  }
}

repeatTP <- filter(PA_TP,uniqueID %in% keepTP)
repeatTN <- filter(PA_TN,uniqueID %in% keepTN)

n = dim(repeatTP)[2]
mutate(repeatTP,mean = rowMeans(repeatTP[,1:(n-1)])) -> repeatTP

m = dim(repeatTN)[2]
mutate(repeatTN,mean = rowMeans(repeatTN[,1:(m-1)])) -> repeatTN

## fill in names
for (i in 1:dim(repeatTN)[1]) {
  for (j in 1:dim(cmpdN)[1]) {
    if (repeatTN$uniqueID[i] == cmpdN$UniqueID[j]) {
      repeatTN$EnglisName[i] <- cmpdN$`English Name`[j]
    }
  }
}


for (i in 1:dim(repeatTP)[1]) {
  for (j in 1:dim(cmpdP)[1]) {
    if (repeatTP$uniqueID[i] == cmpdP$UniqueID[j]) {
      repeatTP$EnglisName[i] <- cmpdP$`English Name`[j]
    }
  }
}

## compare
rpTP <- select(repeatTP,uniqueID:EnglisName)
names(rpTP) <- c("uniqueIDP","meanP","EnglisName")

rpTN <- select(repeatTN,uniqueID:EnglisName)
names(rpTN) <- c("uniqueIDN","meanN","EnglisName")

rp <- merge(rpTN,rpTP,by = "EnglisName")

rpTP <- filter(rp, meanP > meanN)
filter(repeatTP,!(uniqueID %in% rpTP$uniqueIDP)) -> repeatTP

rpTN <- filter(rp, meanN > meanP)
filter(repeatTN,!(uniqueID %in% rpTN$uniqueIDN)) -> repeatTN

TissueP <- filter(TissueP, UniqueID %in% cmpdP$UniqueID)  %>% 
  filter(.,!(UniqueID %in% repeatTP$uniqueID))

TissueN <- filter(TissueN, UniqueID %in% cmpdN$UniqueID) %>% 
  filter(.,!(UniqueID %in% repeatTN$uniqueID))

all(rpTP$uniqueIDP %in% TissueP$UniqueID)
all(rpTN$uniqueIDN %in% TissueN$UniqueID)

####################
#####         ######
##### plasma  ######
#####         ######
####################
## pos
path = "inputFile_V2/Plasma_pos_IS_adj.xlsx"
PlasmaP <- read_excel(path = path,sheet = 1,col_names = T)
PlasmaP <- select(PlasmaP,-(Compounds:QC_injectorder))

## neg
path1 = "inputFile_V2/Plasma_neg_IS_adj.xlsx"
PlasmaN <- read_excel(path = path1,sheet = 1,col_names = T)
PlasmaN <- select(PlasmaN,-(Compounds:QC_inject_order))

rm(path,path1,path2)
load("plasma_medinfo.RData")

PlasmaP1 <- select(PlasmaP,-(QC_1:QC_7)) %>% 
  select(.,UniqueID,contains(medInfoPlasma$sampleID))

PlasmaN1 <- select(PlasmaN,-(QC_1:QC_7)) %>%
  select(.,UniqueID,contains(medInfoPlasma$sampleID))

## new peak area normalization matrix 
## pos
n = dim(PlasmaP1)[2]
PA_PP <- matrix(0,nrow = dim(PlasmaP1)[1],ncol = (n-1),
                dimnames = list(PlasmaP1$UniqueID,names(PlasmaP1)[2:n])) %>% 
  as.data.frame

for (j in 1:(n-1)) {
  for (i in 1:dim(PlasmaP1)[1]) {
    PA_PP[i,j] <- PlasmaP1[i,(j+1)]/sum(PlasmaP1[,(j+1)])
  }
}

PA_PP$uniqueID <- rownames(PA_PP)

## neg
rm(i,j,n)
n = dim(PlasmaN1)[2]
PA_PN <- matrix(0,nrow = dim(PlasmaN1)[1], ncol = (n-1),
                dimnames = list(PlasmaN1$UniqueID,names(PlasmaN1)[2:n])) %>% 
  as.data.frame

for (j in 1:(n-1)) {
  for (i in 1:dim(PlasmaN1)[1]) {
    PA_PN[i,j] <- PlasmaN1[i,(j+1)]/sum(PlasmaN1[,(j+1)])
  }
}

PA_PN$uniqueID <- rownames(PA_PN)

## pick repeated chemicals
rm(i,j,n)
cmpdP <- filter(cmpd,Tissue == "Plasma",adduct_type == "pos")
cmpdN <- filter(cmpd,Tissue == "Plasma",adduct_type == "neg")

keepPP <- NA
keepPN <- NA

for (i in 1:dim(cmpdP)[1]) {
  for (j in 1:dim(cmpdN)[1]) {
    if (cmpdP$`English Name`[i] == cmpdN$`English Name`[j]) {
      keepPP <- append(keepPP,cmpdP$UniqueID[i],0)
      keepPN <- append(keepPN,cmpdN$UniqueID[j],0)
    }
  }
}

repeatPP <- filter(PA_PP,uniqueID %in% keepPP)
repeatPN <- filter(PA_PN,uniqueID %in% keepPN)

n = dim(repeatPP)[2]
mutate(repeatPP,mean = rowMeans(repeatPP[,1:(n-1)])) -> repeatPP

m = dim(repeatPN)[2]
mutate(repeatPN,mean = rowMeans(repeatPN[,1:(m-1)])) -> repeatPN

## fill in names
for (i in 1:dim(repeatPN)[1]) {
  for (j in 1:dim(cmpdN)[1]) {
    if (repeatPN$uniqueID[i] == cmpdN$UniqueID[j]) {
      repeatPN$EnglisName[i] <- cmpdN$`English Name`[j]
    }
  }
}

for (i in 1:dim(repeatPP)[1]) {
  for (j in 1:dim(cmpdP)[1]) {
    if (repeatPP$uniqueID[i] == cmpdP$UniqueID[j]) {
      repeatPP$EnglisName[i] <- cmpdP$`English Name`[j]
    }
  }
}


## compare
rpPP <- select(repeatPP,uniqueID:EnglisName)
names(rpPP) <- c("uniqueIDP","meanP","EnglisName")

rpPN <- select(repeatPN,uniqueID:EnglisName)
names(rpPN) <- c("uniqueIDN","meanN","EnglisName")

rp <- merge(rpPN,rpPP,by = "EnglisName")

rpPP <- filter(rp, meanP > meanN)
filter(repeatPP,!(uniqueID %in% rpPP$uniqueIDP)) -> repeatPP

rpPN <- filter(rp, meanN > meanP)
filter(repeatPN,!(uniqueID %in% rpPN$uniqueIDN)) -> repeatPN

PlasmaP <- filter(PlasmaP, UniqueID %in% cmpdP$UniqueID)  %>% 
  filter(.,!(UniqueID %in% repeatPP$uniqueID))

PlasmaN <- filter(PlasmaN, UniqueID %in% cmpdN$UniqueID) %>% 
  filter(.,!(UniqueID %in% repeatPN$uniqueID))
  
####################
#####         ######
##### URINE   ######
#####         ######
####################
## pos
path = "inputFile_V2/Urine_pos_IS_creatinine_adj.xlsx"
UrineP <- read_excel(path = path,sheet = 1,col_names = T)
UrineP <- select(UrineP,-(Compounds:rtmed))

## neg
path1 <- "inputFile_V2/Urine_neg_IS_creatinine_adj.xlsx"
UrineN <- read_excel(path = path1,sheet = 1, col_names = T)
UrineN <- select(UrineN,-(Compounds:rtmed))

rm(path,path1,path2)
load("Urine_medinfo.RData")
UrineP <- select(UrineP,-(qc.1:qc.7)) %>% 
  select(.,UnqiueID,contains(medInfoUrine$sampleID)) 

UrineN <- select(UrineN,-(qc.1:qc.7)) %>% 
  select(.,UniqueID,contains(medInfoUrine$sampleID))

## new peak area normalization matrix 
## pos
n = dim(UrineP)[2]
PA_UP <- matrix(0,nrow = dim(UrineP)[1],ncol = (n-1),
                dimnames = list(UrineP$UnqiueID,names(UrineP)[2:n])) %>% 
  as.data.frame

for (j in 1:(n-1)) {
  for (i in 1:dim(UrineP)[1]) {
    PA_UP[i,j] <- UrineP[i,(j+1)]/sum(UrineP[,(j+1)])
  }
}

PA_UP$uniqueID <- rownames(PA_UP)

## neg
rm(i,j,n)
n = dim(UrineN)[2]
PA_UN <- matrix(0,nrow = dim(UrineN)[1], ncol = (n-1),
                dimnames = list(UrineN$UniqueID,names(UrineN)[2:n])) %>% 
  as.data.frame

for (j in 1:(n-1)) {
  for (i in 1:dim(UrineN)[1]) {
    PA_UN[i,j] <- UrineN[i,(j+1)]/sum(UrineN[,(j+1)])
  }
}

PA_UN$uniqueID <- rownames(PA_UN)

## pick repeated chemicals
rm(i,j,n)
cmpdP <- filter(cmpd,Tissue == "Urine",adduct_type == "pos")
cmpdN <- filter(cmpd,Tissue == "Urine",adduct_type == "neg")

keepUP <- NA
keepUN <- NA

for (i in 1:dim(cmpdP)[1]) {
  for (j in 1:dim(cmpdN)[1]) {
    if (cmpdP$`English Name`[i] == cmpdN$`English Name`[j]) {
      keepUP <- append(keepUP,cmpdP$UniqueID[i],0)
      keepUN <- append(keepUN,cmpdN$UniqueID[j],0)
    }
  }
}

repeatUP <- filter(PA_UP,uniqueID %in% keepUP)
repeatUN <- filter(PA_UN,uniqueID %in% keepUN)

n = dim(repeatUP)[2]
mutate(repeatUP,mean = rowMeans(repeatUP[,1:(n-1)])) -> repeatUP

m = dim(repeatUN)[2]
mutate(repeatUN,mean = rowMeans(repeatUN[,1:(m-1)])) -> repeatUN

## fill in names
for (i in 1:dim(repeatUN)[1]) {
  for (j in 1:dim(cmpdN)[1]) {
    if (repeatUN$uniqueID[i] == cmpdN$UniqueID[j]) {
      repeatUN$EnglisName[i] <- cmpdN$`English Name`[j]
    }
  }
}

for (i in 1:dim(repeatUP)[1]) {
  for (j in 1:dim(cmpdP)[1]) {
    if (repeatUP$uniqueID[i] == cmpdP$UniqueID[j]) {
      repeatUP$EnglisName[i] <- cmpdP$`English Name`[j]
    }
  }
}

## compare
rpUP <- select(repeatUP,uniqueID:EnglisName)
names(rpUP) <- c("uniqueIDP","meanP","EnglisName")

rpUN <- select(repeatUN,uniqueID:EnglisName)
names(rpUN) <- c("uniqueIDN","meanN","EnglisName")

rp <- merge(rpUN,rpUP,by = "EnglisName")

rpUP <- filter(rp, meanP > meanN)
filter(repeatUP,!(uniqueID %in% rpUP$uniqueIDP)) -> repeatUP

rpUN <- filter(rp, meanN > meanP)
filter(repeatUN,!(uniqueID %in% rpUN$uniqueIDN)) -> repeatUN

UrineP1 <- filter(UrineP, UnqiueID %in% cmpdP$UniqueID)  %>% 
  filter(.,!(UnqiueID %in% repeatUP$uniqueID))

UrineN1 <- filter(UrineN, UniqueID %in% cmpdN$UniqueID) %>% 
  filter(.,!(UniqueID %in% repeatUN$uniqueID))

all(rpUP$uniqueIDP %in% UrineP1$UnqiueID)
all(rpUN$uniqueIDN %in% UrineN1$UniqueID)


