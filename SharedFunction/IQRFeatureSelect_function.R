### Input -- sample in row, features in column
## first two columns: smapleID & TCP
## column-based outliers(variables-based)

featureSelect <- function(mydata) {
  n = ncol(mydata)
  m = nrow(mydata)
  
  for (i in 3:n) {
    b = table(mydata[1:m,i]) %>% as.data.frame
    a = which(b[,1] == "outliers")
    c = b[a,2]
    
    if (length(c) == 0) {
      mydata[(m+1),i] <- "keep"
      }else if (c/m > 0.7) {
        mydata[(m+1),i] <- "unkeep"
        } else if (c/m < 0.7) {
          mydata[(m+1),i] <- "keep"
          }
    }
  
  return(mydata)	
}

save(featureSelect, file = "featureSelect_Func.RData")
