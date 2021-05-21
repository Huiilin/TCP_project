### Input -- sample in row, features in column
# first two columns: smapleID & TCP

outliers2NA <- function(mydata) {
	m = nrow(mydata)
	n = ncol(mydata)
	a = which(names(mydata) == "TCP")

	for (i in (a+1):n) {
		for (j in 1:(m-1)) {
			if (is.na(mydata[j,i]) == "TRUE") {next};
			if (mydata[j,i] == "outliers") mydata[j,i] <- NA
		}
	}

	return(mydata[(1:(m-1)),])
}

save(outliers2NA, file = "outliers2NAFunc.RData")
