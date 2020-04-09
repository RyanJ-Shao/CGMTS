# The MIT License (MIT)
# Copyright (c) 2020 UCAS

#Use manual format
prepro <- function(inputdir="", outputdir="", outlierdet = TRUE, interval = 15, imputation = FALSE, 
                   immethod = "linear", maxgap = 60, compeleteday = TRUE, removeday = FALSE){
	fileNames = list.files(inputdir)
	for(f in fileNames){
	  print(paste("processing file:", f))
	  cgmts <- read.csv(paste(inputdir, "/", f, sep = ''))
	  colnm <- colnames(cgmts)
	  if(colnm[1] != "timestmap" || colnm[2] != "sglucose" || colnm[3] != "bglucose"){
	    stop(paste("The format fo file '",f ,"' is incorrect and cannot be read.",sep = ""))
	  }
	  cgmts  = qcfun(cgmts, outlierdet, interval, imputation,immethod, maxgap, compeleteday,removeday)
	  write.csv(cgmts, paste(outputdir,"/",f,sep=''),row.names = FALSE)
	  }
}
