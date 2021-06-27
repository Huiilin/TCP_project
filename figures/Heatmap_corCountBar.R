rm(list = ls())
options(scipen = 200)

library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

path = "15_heatmap_mofa_imputedData/combine_LF1weight_cmpd_NoSynthetic_20210602.xlsx"
weight <- read_excel(path,sheet = 1,col_names = T)
impute <- read_excel("15_heatmap_mofa_imputedData/imputedData_withLF1weight_1.xlsx",
                     sheet = 1)
                     
plasma <- filter(weight, view == "Plasma") %>% 
  distinct(featureID,.keep_all = TRUE) %>% arrange(by_group = desc(absWeight))

plasma26 <- filter(plasma,absWeight > 0.4) %>% 
  dplyr::select("featureID","Annotation","absWeight")

imputePlasma <- filter(impute,view == "Plasma")

## combine impute data & weight data
plasma26impute <- left_join(plasma26,imputePlasma,by = "featureID")
plasma26impute <- dplyr::select(plasma26impute,-(Annotation.x:upDown_1vs2)) %>%
  as.data.frame()
rownames(plasma26impute) <- plasma26impute$featureID

plasmaHeatmap <- plasma26impute[,-1] %>% as.matrix()

## read in sample info 
load("1_clinical_info/sampleInfoAll.RData")
names(sampleInfo)
Group <- sampleInfo[,c("病案号","血小板减少")]
names(Group) <- c("medicalID","TCP")
which(Group$medicalID  == "741044") 
Group$medicalID[41] <- '741044_1'

Group <- filter(Group,medicalID %in% colnames(plasmaHeatmap))
Group$TCP <- as.factor(Group$TCP)
plasmaHeatmap1 <- plasmaHeatmap[,Group$medicalID]

f1 <- colorRamp2(seq(min(plasmaHeatmap1), max(plasmaHeatmap1), length=3), c("blue","#EEEEEE", "red"))

ha2 <- HeatmapAnnotation(TCP = Group$TCP,
                         col = list(TCP = c("0" = "#2ca25f","1" = "#fec44f","2" = "#de2d26")))

p1 <-Heatmap(plasmaHeatmap1,cluster_rows = T,
             column_split = Group$TCP,
             col = f1,name = "value",color_space = "LAB",
             show_column_names = T,
             top_annotation = ha2)
p1

## add barplot
## read in all counts data
counts <- read.csv("14_corCont_all/count_cor_all.csv",header = T)
plasma26counts <- filter(counts,featureID %in% plasma26$featureID) %>% .[,-2] 

## re-order count data as the same order of heatmap rows
a = row_order(p1)
aa = rownames(plasmaHeatmap1)[a] %>% as.data.frame()
names(aa) <- "featureID"

plasma26counts1 <- left_join(aa,plasma26counts,by = "featureID")
rownames(plasma26counts1) <- plasma26counts1$featureID
plasma26counts1 <- plasma26counts1[,-1] %>% as.matrix()

mycolor <- c("1" = "#2c7fb8","2" = "#feb24c","3" = "#f03b20","4" = "#756bb1","5" = "#2ca25f")
  
## log2 counts
plasma26counts1 <- log2(plasma26counts1)

ha3 <- rowAnnotation(cor = anno_barplot(plasma26counts1,
                                        gp = gpar(fill = mycolor),
                                        width = unit(25,"cm"),
                                        border = F,
                                        axis = T,
                                        bar_width = 0.8))
                                        
p2 <- Heatmap(plasmaHeatmap1,cluster_rows = T,
        column_split = Group$TCP,
        col = f1,name = "value",color_space = "LAB",
        show_column_names = F,
        show_row_names = T,
        top_annotation = ha2,
        right_annotation = ha3)
p2

## anno_barplot of mofa absweight
ha4color <- c("1" = '#e7e1ef')
ha4 <- rowAnnotation(weight = anno_barplot(plasma26_1[,3],
                                           gp = gpar(fill = ha4color),
                                           axis_param = list(direction = "reverse"),
                                           width = unit(25,"cm"),
                                           border = F,
                                           axis = T,
                                           bar_width = 0.8,
                                           height = unit(6,"cm")))
p3 <- draw(ha4 + p2)
p3

## add legend
lgd <- Legend(at = colnames(plasma26counts1)[1:5],
              legend_gp = gpar(fill = mycolor),
              labels = c("Plasma","Urine","Tissue","RNA","Methylation"),
              type = "grid",
              title = "Omics Groups")

lgd2 <- Legend(at = colnames(plasma26_1)[3],legend_gp = gpar(fill = ha4color),
               title = "absolute weight of MOFA weights",
               type = "grid")

pd <- packLegend(lgd,lgd2)

## with legend
p4 <- draw(p2,annotation_legend_list = pd)
p4

## save 
library(export)
graph2ppt(p3, "15_heatmap_plasma26_corFigure/heatmap_cor_v1.pptx",append = TRUE,
          width = 45,height = 15)

graph2ppt(p4, "15_heatmap_plasma26_corFigure/heatmap_cor_v1.pptx",append = TRUE,
          width = 30,height = 15)
  
