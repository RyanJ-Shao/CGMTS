# The MIT License (MIT)
# Copyright (c) 2020 UCAS
library(dplyr)
library(imputeTS)
library(forecast)
library(TSA)
qcfun<- function(cgmts, outlierdet = TRUE, interval = 15, imputation = FALSE, immethod = "linear", 
                 maxgap = 60, compeleteday = TRUE, removeday = FALSE, units = 'mmol/L'){
  vectimestamp <- as.vector(cgmts$timestamp)
  vectimestamp <- unlist(strsplit(vectimestamp,split=" "))
  maxtimestamp <- matrix(vectimestamp,ncol=2,byrow=T)[,1]
  cgmts <- cgmts %>% mutate(timedate = maxtimestamp)
  coldate <- unique(cgmts$timedate)
  #remove first day and last day
  fday <- coldate[1]
  lday <- coldate[-1]
  #some data points missed, the time still be recored
  if(length(cgmts[cgmts$timedate ==fday,]$timedate) <freq){
    cgmts[cgmts$timedate ==fday,]$timedate <- NA
  }
  if(length(cgmts[cgmts$timedate ==lday,]$timedate) <freq){
    cgmts[cgmts$timedate ==lday,]$timedate <- NA
  }
  cgmts <- cgmts[!is.na(cgmts$timedate),]
  if(units != "mmol/L"){
    cgmts <- mutate(cgmts, sglucose = round(sglucose/18,2))
  }
  freq = 1440/interval
  if(compeleteday){
    for (d in coldate){
      if (length(cgmts[cgmts$timedate == d,]$timedate) < freq){
        cgmts[cgmts$timedate ==d,]$timedate <- NA
      }
    }
    cgmts <- cgmts[!is.na(cgmts$timedate),]
  }else if(imputation){
    if(immethod == "linear"){
      print("linear imputation")
      cgmts <- cgmts %>% mutate(imglucose = na.interpolation(cgmts$sglucose,maxgap = as.integer(maxgap/interval)))
    }else if(immethod == "seadec"){
      print("SEADEC imputation")
      tstype <- ts(cgmts$sglucose,frequency = freq)
      tstype <- na.seadec(tstype,algorithm = "interpolation",maxgap = as.integer(maxgap/interval))
      cgmts <- cgmts %>% mutate(imglucose = tstype)
    }else if(immethod == "arima"){
      print("ARIMA imputation")
      cgmts <- cgmts %>% mutate(imglucose = na.kalman(cgmts$sglucose,model = "auto.arima",maxgap = as.integer(maxgap/interval)))
      
    }
    if(removeday == TRUE){
      is.na.rle <- rle(is.na(cgmts$sglucose))
      is.na.rle$values <- is.na.rle$values & is.na.rle$lengths > as.integer(maxgap/interval)
      reov = cgmts[inverse.rle(is.na.rle), ]
      reovdate = unique(reov$timedate)
      for (d in reovdate){
        cgmts[cgmts$timedate ==d,]$timedate <- NA
      }
      cgmts <- cgmts[!is.na(cgmts$timedate),]
    }
  }else{
    if(removeday == TRUE){
      is.na.rle <- rle(is.na(cgmts$sglucose))
      is.na.rle$values <- is.na.rle$values & is.na.rle$lengths > as.integer(maxgap/interval)
      reov = cgmts[inverse.rle(is.na.rle), ]
      reovdate = unique(reov$timedate)
      for (d in reovdate){
        cgmts[cgmts$timedate ==d,]$timedate <- NA
      }
      cgmts <- cgmts[!is.na(cgmts$timedate),]
  }
  
  }
  if(outlierdet == TRUE){
    unidate = unique(cgmts$timedate)
    for (d in unidate){
      udcgm <- cgmts[cgmts$timedate ==d,]
      if(!any(is.na(udcgm$sglucose))){
        udcgmglucose <- udcgm$sglucose
        udcgmts <- ts(udcgmglucose, frequency = freq)
        tsmodel <- auto.arima(udcgmts)
        ao <- detectAO(tsmodel)
        io <- detectIO(tsmodel)
        aodict <- c(ao$lambda2)
        names(aodict) <- ao$ind
        
        iodict <- c(io$lambda1)
        names(iodict) <- io$ind
        indinc <- intersect(ao$ind, io$ind)
        #remove ao or io index according to theri lambda
        for(ic in indinc){
          if(aodict[[as.character(ic)]] < iodict[[as.character(ic)]]){
            removeind = which(c(ic) %in% names(aodict) )
            aodict <- aodict[-removeind]
          }else{
            removeind <- which(c(ic) %in% names(iodict) )
            iodict <- iodict[-removeind]
          }
        }
        
        for(i in seq_along(aodict)){
          ind <- as.numeric(names(aodict)[i])
          udcgm[c(ind),]$outliers <- "AO"
        }
        for(i in seq_along(iodict)){
          ind <- as.numeric(names(iodict)[i])
          udcgm[c(ind),]$outliers <- "IO"
        }
        cgmts[cgmts$timedate ==d,]$outliers <- udcgm$outliers
      }
    }
    
    return(cgmts)
  }
}