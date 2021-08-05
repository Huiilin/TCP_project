rm(list = ls())
options(scipen = 200)

load("gene_FPKM_54_KNNimputed.RData")
load("Trans_medinfo.RData")

## set clinical sample
keep <- colnames(imputeDat2)
RNA <- medInfoTrans[keep,]

## adjust for covariates
design <- model.matrix(~0+condition+Sex+ageGroup+HeightGroup+WeightGroup
                       +Staging+BMIGroup,data = RNA)
nonEstimable(design)
RNA_Adjust <- removeBatchEffect(imputeDat2,batch = RNA$HistologyType,
                              design = design)

design <- model.matrix(~0+condition+ageGroup+HeightGroup+WeightGroup
                       +Staging+BMIGroup,data = RNA)
RNA_Adjust <- removeBatchEffect(RNA_Adjust,batch = RNA$Sex,design = design)

design <- model.matrix(~0+condition+HeightGroup+WeightGroup+Staging+BMIGroup,
                       data =  RNA)
RNA_Adjust <- removeBatchEffect(RNA_Adjust,batch = RNA$ageGroup,design = design)

design <- model.matrix(~0+condition+WeightGroup+Staging+BMIGroup,data = RNA)
RNA_Adjust <- removeBatchEffect(RNA_Adjust,batch = RNA$HeightGroup,design = design)

design <- model.matrix(~0+condition+Staging+BMIGroup,data = RNA)
RNA_Adjust <- removeBatchEffect(RNA_Adjust,batch = RNA$WeightGroup,design = design)

design <- model.matrix(~0+condition+BMIGroup,data = RNA)
RNA_Adjust <- removeBatchEffect(RNA_Adjust,batch = RNA$Staging,design = design)

design <- model.matrix(~0+condition,data = RNA)
RNA_Adjust <- removeBatchEffect(RNA_Adjust,batch = RNA$BMIGroup,design = design)

## DE
contrMatrix <- makeContrasts(TCP0vs1 = conditionTCP1-conditionNoTCP,
                             TCP0vs2 = conditionTCP2-conditionNoTCP,
                             TCP1vs2 = conditionTCP2-conditionTCP1,
                             levels = colnames(design))

fit <- lmFit(RNA_Adjust,design = design)
fit <- contrasts.fit(fit = fit,contrasts = contrMatrix)
efit <- eBayes(fit)

## get results
res10 <- topTable(efit,coef = 1, number = nrow(RNA_Adjust))
diffMetabos10 <- res10[which(res10$P.Value < 0.05),]
ID10 <- rownames(diffMetabos10)

res20 <- topTable(efit,coef = 2, number = nrow(RNA_Adjust))
diffMetabos20 <- res20[which(res20$P.Value < 0.05),]
ID20 <- rownames(diffMetabos20)

res21 <- topTable(efit,coef = 3, number = nrow(RNA_Adjust))
diffMetabos21 <- res21[which(res21$P.Value < 0.05),]
ID21 <- rownames(diffMetabos21)

IDAll <- union(ID10,ID20) %>% union(.,ID21)

RNA_DE <- RNA_Adjust[IDAll,] %>% as.data.frame

###
save(RNA_DE,RNA,file = "gene_FPKM_54_DE.RData")
