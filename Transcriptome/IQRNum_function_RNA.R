# input -- sample in row, feature in column
# first two columns: smapleID & TCP
## variables based IQR
IQRNum <- function(mydata) {
	cols <- str_subset(colnames(mydata), pattern = "ENSG")

	n = nrow(mydata)
	mydata[(n+1),cols] = apply(mydata[,cols],2,IQR,na.rm =T)
	rownames(mydata)[(n+1)] <- "IQR"
  
	return(mydata)
}

save(IQRNum, file = "IQRNumFunc.RData")
