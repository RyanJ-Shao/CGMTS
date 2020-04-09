# The MIT License (MIT)
# Copyright (c) 2020 UCAS
library(dplyr)
library(plotly)

cgmplot <- function(inputdir, outputdir, useig= FALSE, markoutliers= TRUE, interval = 15,
                    diffnum = 1, seadiff = TRUE, html = TRUE){
  fileNames = list.files(inputdir)
  for(f in fileNames){
    fname = unlist(strsplit(f, split = "\\."))[1]
    print(paste("processing file:", f))
    cgmtsall = read.csv(paste(inputdir, "/", f, sep = ''),stringsAsFactors= FALSE)
    vectimestamp <- as.vector(cgmtsall$timestamp)
    vectimestamp <- unlist(strsplit(vectimestamp,split=" "))
    maxtimestamp <- matrix(vectimestamp,ncol=2,byrow=T)[,1]
    cgmtsall <- cgmtsall %>% mutate(timedate = maxtimestamp)
    #print(head(cgmtsall))
    print("plottinig ACF")
    cgmACF(cgmtsall, fname, outputdir, useig, diffnum, seadiff, interval)
    print("plottinig PACF")
    cgmPACF(cgmtsall, fname, outputdir, useig, diffnum, seadiff, interval)
    print("plottinig 3d")
    cgm3d(cgmtsall, fname, outputdir, useig)
    print("plottinig decom")
    cgmdecom(cgmtsall, fname, outputdir, useig,interval, html = html)
    print("plottinig trace")
    cgmtrace(cgmtsall, fname, outputdir,useig, markoutliers, html = html)
  }
  
  
}


cgmACF <- function(cgmtsall, fname, outputdir, useig = TRUE, diffnum = 1, seadiff = TRUE, interval = 15){
  glucosets = NULL
  if(useig){
    glucosets = cgmtsall$imglucose
  }else{
    glucosets = cgmtsall$sglucose
  }
  glucosets = ts(glucosets, frequency = 1440/15)
  pdf(paste(outputdir,"/", fname,"_","acf", ".pdf",sep = "")) 
  if(seadiff){
    acf(diff(glucosets, differences = diffnum), lag = 1440/15)
  }else{
    acf(diff(glucosets, differences = diffnum))
  }
  
  dev.off()
}


cgmPACF <- function(cgmtsall, fname, outputdir, useig = TRUE, diffnum = 1, seadiff = TRUE, interval = 15){
  glucosets = NULL
  if(useig){
    glucosets = cgmtsall$imglucose
  }else{
    glucosets = cgmtsall$sglucose
  }
  glucosets = ts(glucosets, frequency = 1440/interval)
  pdf(paste(outputdir,"/", fname,"_","pacf", ".pdf",sep = "")) 
  if(seadiff){
    pacf(diff(glucosets, differences = diffnum), lag = 1440/interval)
  }else{
    pacf(diff(glucosets, differences = diffnum))
  }
  dev.off()
}


cgm3d <- function(cgmtsall, fname, outputdir, useig = TRUE){
  #x = col, y row
  x = c('00:00:00','00:15:00','00:30:00','00:45:00','01:00:00','01:15:00','01:30:00',
        '01:45:00','02:00:00','02:15:00','02:30:00','02:45:00','03:00:00','03:15:00',
        '03:30:00','03:45:00','04:00:00','04:15:00','04:30:00','04:45:00','05:00:00',
        '05:15:00','05:30:00','05:45:00','06:00:00','06:15:00','06:30:00','06:45:00',
        '07:00:00','07:15:00','07:30:00','07:45:00','08:00:00','08:15:00','08:30:00',
        '08:45:00','09:00:00','09:15:00','09:30:00','09:45:00','10:00:00','10:15:00',
        '10:30:00','10:45:00','11:00:00','11:15:00','11:30:00','11:45:00','12:00:00',
        '12:15:00','12:30:00','12:45:00','13:00:00','13:15:00','13:30:00','13:45:00',
        '14:00:00','14:15:00','14:30:00','14:45:00','15:00:00','15:15:00','15:30:00',
        '15:45:00','16:00:00','16:15:00','16:30:00','16:45:00','17:00:00','17:15:00',
        '17:30:00','17:45:00','18:00:00','18:15:00','18:30:00','18:45:00','19:00:00',
        '19:15:00','19:30:00','19:45:00','20:00:00','20:15:00','20:30:00','20:45:00',
        '21:00:00','21:15:00','21:30:00','21:45:00','22:00:00','22:15:00','22:30:00',
        '22:45:00','23:00:00','23:15:00','23:30:00','23:45:00')
  if(useig){
    z = matrix(cgmtsall$imglucose,ncol=96,byrow=T)
  }else{
    z = matrix(cgmtsall$sglucose,ncol=96,byrow=T)
  }
  fig <- plot_ly(x = x, z = z)
  fig <- fig %>% add_surface()
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = list(
        title = "Time",
        dtick = 10, 
        tick0 = 0
        #tickmode = "array",
        #type = "date",
        #tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Days"
      ),
      zaxis = list(title = "Glucose Values")
    )
    )
  
  htmlwidgets::saveWidget(fig, paste(outputdir,"/", fname,"_","3d", ".html",sep = ""))
}


cgmdecom <- function(cgmtsall, fname, outputdir, useig = TRUE,interval = 15, html = FALSE){
  freq = 1440/interval
  uniday = unique(cgmtsall$timedate)
  #remove uncompelete day
  for(d in uniday){
    if(useig){
      if(any(is.na(cgmtsall[cgmtsall$timedate ==d,]$imglucose))){
        cgmtsall[cgmtsall$timedate ==d,]$timedate <- NA
      }
    }else{
      if(any(is.na(cgmtsall[cgmtsall$timedate ==d,]$sglucose))){
        cgmtsall[cgmtsall$timedate ==d,]$timedate <- NA
      }
    }
  }
  cgmtsall <- cgmtsall[!is.na(cgmtsall$timedate),]
  uniday = unique(cgmtsall$timedate)
  if(useig){
    gts <- ts(cgmtsall$imglucose, frequency = freq)
  }else{
    gts <- ts(cgmtsall$sglucose, frequency = freq)
  }
  stlgts <- stl(gts,s.window = "periodic")
  seasonal <- stlgts$time.series[,1]
  trend <- stlgts$time.series[,2]
  remainder <- stlgts$time.series[,3]
  #plot seasonal
  seafig <- plot_ly(        
    type = "scatter",
    x = c(cgmtsall$timestamp),
    y = c(seasonal),
    mode = "lines"
  )
  seafig <- seafig %>%
    layout(
      xaxis = list(
        title = "Time",
        dtick = 10, 
        tick0 = 0,
        tickmode = "array",
        type = "date",
        tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Glucose seansonal component"
      )
    )
  #plot trend
  trfig <- plot_ly(        
    type = "scatter",
    x = c(cgmtsall$timestamp),
    y = c(trend),
    mode = "lines"
  )
  trfig <- trfig %>%
    layout(
      xaxis = list(
        title = "Time",
        dtick = 10, 
        tick0 = 0,
        tickmode = "array",
        type = "date",
        tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Glucose trend component"
      )
    )

  #plot trend
  refig <- plot_ly(        
    type = "scatter",
    x = c(cgmtsall$timestamp),
    y = c(remainder),
    mode = "lines"
  )
  refig <- refig %>%
    layout(
      xaxis = list(
        title = "Time",
        dtick = 10, 
        tick0 = 0,
        tickmode = "array",
        type = "date",
        tickformat = "%H:%M:%S<br>%Y-%B-%d"
      ),
      yaxis = list(
        title = "Glucose remainder component"
      )
    )

  if(html){
    htmlwidgets::saveWidget(seafig, paste(outputdir,"/", fname,"_","seasonal", ".html",sep = ""))
    htmlwidgets::saveWidget(trfig, paste(outputdir,"/", fname,"_","trend", ".html",sep = ""))
    htmlwidgets::saveWidget(refig, paste(outputdir,"/", fname,"_","remainder", ".html",sep = ""))
  }else{
    oldworkdir = getwd()
    setwd(outputdir)
    orca(seafig, paste(fname,"_","seasonal", ".pdf",sep = ""))
    orca(trfig, paste( fname,"_","trend", ".pdf",sep = ""))
    orca(refig, paste(fname,"_","remainder", ".pdf",sep = ""))
    setwd(outputdir)
  }
  
  #orca(seafig, paste(outputdir,"/", fname,"_","seasonal", ".pdf",sep = ""))
  #orca(trfig, paste(outputdir,"/", fname,"_","trend", ".pdf",sep = ""))
  #orca(refig, paste(outputdir,"/", fname,"_","remainder", ".pdf",sep = ""))

  #return(seafig) 
}




#fname = ryan
#cgmtrace(cgmts, "ryan", "Desktop/cgm_software/CGMTS/plotoutput")
cgmtrace <- function(cgmtsall, fname, outputdir,useig = TRUE, markoutliers = TRUE, html = FALSE){

  cgmdate <- unique(cgmtsall$timedate)
  for (d in cgmdate){
    cgmts <-filter(cgmtsall, timedate == d)
    vectimestamp <- as.vector(cgmts$timestamp)
    vectimestamp <- unlist(strsplit(vectimestamp,split=" "))
    hms <- matrix(vectimestamp,ncol=2,byrow=T)[,2]
    gts <- c()
    gy <-hms
    if(useig){
      gts <- cgmts$imglucose
    }else{
      gts <- cgmts$sglucose
    }

    fig <- plot_ly(        
      type = "scatter",
      x = c(gy),
      y = c(gts),
      name = 'Glucose Trace',
      mode = "lines"
      )

    if(markoutliers){
      out <- filter(cgmts, outliers == "IO" |  outliers == "AO" )
      for(i in 1:nrow(out)){
        y=0
        x = ""
        if(useig){
          y = out[i,]$imglucose
        }else{
          y = out[i,]$sglucose
        }
        fig <- fig %>%
          add_trace(
            type = "scatter",
            x = unlist(strsplit(out[i,]$timestamp, split = " "))[2],
            y = y,
            name = out[i,]$outliers,
            mode = "markers"
          ) 
      }
    }

    
    fig <- fig %>%
      layout(
        xaxis = list(
          title = "Time",
          dtick = 10, 
          tick0 = 0
        ),
        yaxis = list(
          title = "Glucose value(mmol/L)"
        )
        )
    if(html){
      htmlwidgets::saveWidget(fig, paste(outputdir,"/", fname,"_", d , ".html",sep = ""))
    }else{
      oldworkdir = getwd()
      setwd(outputdir)
      orca(fig, paste(fname, "_", d , ".pdf",sep = ""))
      setwd(oldworkdir)
    }

  }
}