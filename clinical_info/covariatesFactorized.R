options(stringsAsFactors = F)
options(scipen = 200)

library(tidyverse)

## read in files
medInfo <- read.csv(file,na = "NA",header = TRUE)

medInfo <- medInfo[,c(1:5,8:12)]
names(medInfo)
names(medInfo)[4] <- "HistologyType"
names(medInfo)[2] <- "medicalID"

# covariates
## factor HistologyType
medInfo$HistologyType <- factor(medInfo$HistologyType)

## tumor staging 
medInfo$Staging[medInfo$PTNM == "1"] = 1
medInfo$Staging[medInfo$PTNM == "2a"] = 1
medInfo$Staging[medInfo$PTNM == "2b"] = 1
medInfo$Staging[medInfo$PTNM == "2c"] = 1
medInfo$Staging[medInfo$PTNM == "3a"] = 2
medInfo$Staging[medInfo$PTNM == "3b"] = 2
medInfo$Staging[medInfo$PTNM == "3c"] = 2
medInfo$Staging[medInfo$PTNM == "4a"] = 3

medInfo$Staging <- factor(medInfo$Staging)
table(medInfo$PTNM)

## factor sex
medInfo$Sex <- factor(medInfo$Sex)

## factor age, cut at 50,60,70
n = max(medInfo$Age_2019Record) + 1
medInfo$ageGroup <- cut(medInfo$Age_2019Record,breaks = c(1,50,60,70,n),right = F)
table(medInfo$ageGroup)

## factor Height 
a = median(medInfo$Height)
a1 = min(medInfo$Height)
a2 = max(medInfo$Height)
medInfo$HeightGroup <- cut(medInfo$Height, breaks = c(a1,a,a2), right = T,include.lowest = T)
table(medInfo$HeightGroup)

## factor Weight 
a = median(medInfo$Weight)
a1 = min(medInfo$Weight)
a2 = max(medInfo$Weight)
medInfo$WeightGroup <- cut(medInfo$Weight, breaks = c(a1,a,a2), right = T,include.lowest = T)

## factor BMI
a1 = min(medInfo$BMI)
a2 = max(medInfo$BMI)
medInfo$BMIGroup <- cut(medInfo$BMI, breaks = c(a1,18.5,23,a2), right = T,include.lowest = T)

save(medInfo,file = "medInfo_covariatesFactorised.RData")
