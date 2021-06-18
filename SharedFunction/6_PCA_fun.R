## input: de output data,batch effect adjust
## features in rows, sample in column
library(ggord)
library(tidyverse)

metaboPCA <- function(mydata,medinfo) {
  
  a = mydata %>% t(.) %>% as.data.frame
  a$sampleID <- rownames(a)
  
  b = dplyr::select(medinfo,c("TCP","sampleID"))
  
  df = merge(b,a,by = "sampleID")
  df$TCP <- as.factor(df$TCP)
  n = ncol(df)

  pca <- prcomp(df[,3:n])
  
  p = ggord(pca,df$TCP,ploy = FALSE,arrow = 0,ellipse = T,
            vec_ext =0,txt = NULL,ellipse_pro = 0.95)
  p 
  
  return(p)
}
