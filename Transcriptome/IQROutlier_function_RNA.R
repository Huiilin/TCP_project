### Input -- sample in row, features in column
# first two columns: smapleID & TCP
## column-based outliers(variables-based)
outliers <- function(mydata) {
  cols <- str_subset(colnames(mydata), pattern = "ENSG")
  a = which(names(mydata) == "TCP")
  n = ncol(mydata)
  m = nrow(mydata)

  for (i in (a+1):n) {
    Q3 = quantile(mydata[(1:(m-1)),i],0.75,na.rm = TRUE)
    Q1 = quantile(mydata[(1:(m-1)),i],0.25,na.rm = TRUE)

    mydata[(m+1),i] = Q3 + 1.5*mydata[m,i]
    rownames(mydata)[m+1] <- "valueQ3"
    mydata[(m+2),i] = Q1 - 1.5*mydata[m,i]
    rownames(mydata)[m+2] <- "valueQ1"
  }  

  for (i in (a+1):n) {
    for (j in 1:(m-1)) {
      if (is.na(mydata[j,i]) == TRUE) {next}; 
      if (mydata[j,i] < mydata[(m+2),i] | mydata[j,i] > mydata[(m+1),i]) {
        mydata[j,i] <- "outliers"}
    }
  }
  return(mydata)
}


save(outliers, file = "IQROutlierFunc.RData")
