# The MIT License (MIT)
# Copyright (c) 2020 UCAS


#device = 0: manual format
#device = 1: libre freestyle
#device = 2: ipro2
#device = 3: g6
prepro <- function(inputdir="", outputdir="", outlierdet = TRUE, interval = 15, imputation = FALSE,
                   immethod = "linear", maxgap = 60, compeleteday = TRUE, removeday = FALSE, device = 0, transunits = FALSE, removeflday = TRUE){
	fileNames = list.files(inputdir)
	for(f in fileNames){
	  print(paste("processing file:", f))
	  cgmts <- fformat(fpath = paste(inputdir, f, sep = ''), device = device)
	  colnm <- colnames(cgmts)
	  if(colnm[1] != "timestamp" || colnm[2] != "sglucose" || colnm[3] != "bglucose"){
	    stop(paste("The format fo file '",f ,"' is incorrect and cannot be read.",sep = ""))
	  }
	  cgmts  = qcfun(cgmts, outlierdet, interval, imputation,immethod, maxgap, compeleteday,removeday, transunits, removeflday)
	  write.csv(cgmts, paste(outputdir,f,sep=''),row.names = FALSE)
	  }
}

fformat <- function(fpath, device = 0){
  cgmts = ""
  if(device == 0){
    cgmts <- read.csv(fpath)
    return(cgmts)
  }else if(device == 1){
    cgmts <- read.table(fpath, sep = "\t", skip = 2,header = TRUE)
    cgmts <- dplyr::select(cgmts, 2,4)
    names(cgmts) <- c("timestamp", "sglucose")
    cgmts <-  dplyr::mutate(cgmts, timestamp = gsub("/","-", cgmts$timestamp))
    cgmts <-  dplyr::mutate(cgmts, bglucose = NA)
  }else if(device == 2){
    cgmts <- read.table(fpath, sep = "\t", skip = 11,header = TRUE)
    cgmts <- dplyr::select(cgmts, 4,10)
    names(cgmts) <- c("timestamp", "sglucose")
    x <- as.character(lubridate::dmy_hms(cgmts$timestamp))
    cgmts$timestamp <- x
    cgmts <-  dplyr::mutate(cgmts, bglucose = NA)
  }else if(device == 3){
    cgmts <- read.table(fpath, sep = ",", skip = 11)
    cgmts <- dplyr::select(cgmts, 2, 8)
    names(cgmts) <- c("timestamp", "sglucose")
    cgmts <-  dplyr::mutate(cgmts, timestamp = gsub("T"," ", cgmts$timestamp))
    x <- as.character(lubridate::ymd_hms(cgmts$timestamp))
    cgmts$timestamp <- x
    cgmts <-  dplyr::mutate(cgmts, bglucose = NA)
  }
  return(cgmts)
}

