#' General phenology forecast function
#'
#' @param forecastType The type of forecast (randomWalk or logistic)
#' @param forecastLength The length of the forecast into the future in days
#' @param siteName The site name
#' @param URL The URL where the site's phenocam data is located
#' @param lat The latitude of the site in decimals
#' @param long The longitude of the site in decimals
#' @param dataDirectory The file path for the directory to download and store data in
#' @param startDate The start date for the forecast in date form
#' @param endDate The end date for the forecast in date form
#' @param cVals The c values for rescaling
#' @param dVals The d values for rescaling
#' @param GEFS_Files The filenames for the GEFS files
#' @import rjags
#' @import runjags
#' @import coda
#' @import PhenologyBayesModeling
#' @export
phenologyForecast <- function(forecastType,forecastLength=16,siteName,URL,lat,long,dataDirectory,startDate,endDate,season="",cVals,dVals,GEFS_Files){
  nchain=5
  ###Download PhenoCam data and format
  phenoData <- download.phenocam(URL)
  p.old <- phenoData$gcc_mean
  time.old <-  as.Date(phenoData$date)
  days <- seq(as.Date(startDate),(as.Date(endDate)+forecastLength),"day")
  p <- rep(NA,length(days))

  for(i in 1:length(p.old)){
    p[which(days==time.old[i])] <- p.old[i]
  }

  ##Download and format MODIS NDVI data
  mn <- prepareMODIS(startDate=startDate,endDate=endDate,metric="NDVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)
  me <- prepareMODIS(startDate=startDate,endDate=endDate,metric="EVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)

  months <- lubridate::month(days)
  years <- lubridate::year(days)
  data <- list(x_ic=0,tau_ic = 1/(phenoData$g_std[1]**2))

  if(forecastType=="randomWalk"){
    data$p <- rescaleObs(times=days,vals=p,partialStart=TRUE,cVals=cVals,dVals=dVals)
    data$n <- length(data$p)
    data$mn <- rescaleObs(times=days,vals=mn,partialStart=TRUE,cVals=cVals,dVals=dVals)
    data$me <- rescaleObs(times=days,vals=me,partialStart=TRUE,cVals=cVals,dVals=dVals)
    j.model <- randomWalkPhenoModel(data=data,nchain=nchain)
    variableNames <- c("p.PC","p.MN","p.ME","p.proc","x")
  }
  else if(forecastType=="logistic" || forecastType== "logisticCov"){
    dat2 <- data.frame(dates=days,years=years,months=months,p=p,mn=mn,me=me)
    SfsALL <- matrix(nrow=0,ncol=length(dat2$dates))
    if(forecastType=="logisticCov"){
      Tairs <- download_US_WCr_met(start_date=startDate,end_date=endDate)
      for(G in 1:length(GEFS_Files)){
        TairsForecast <- load_GEFS_Forecast(paste("PhenologyForecastData/",GEFS_Files[G],sep=""))
        fileTairs <- c(Tairs,TairsForecast)
        Sfs <- calSf(Tairs=fileTairs,days=as.Date(dat2$dates))
        SfsALL <- rbind(SfsALL,Sfs)
      }
      SfsMeans <- colMeans(SfsALL)
      SfsVar <- apply(SfsALL,MARGIN=2,FUN=var)
      #dat2$SfsPast <- SfsMeans[SfsVar==0]
      dat2$Sf <- SfsMeans
      #SfsFuture <- SfsMeans[(length(SfsMeans)-15):length(SfsMeans)]
      SfsFuture <- SfsMeans[SfsVar!=0]
      #SfsPrec <- 1/SfsVar[(length(SfsVar)-15):length(SfsVar)]
      SfsPrec <- 1/SfsVar[SfsVar!=0]
      print(SfsFuture)
      print(SfsPrec)
    }

    dat2 <- dat2[dat2$months%in%seq(1,6,1),]

    p <- matrix(nrow=181,ncol=0)
    Sf <- matrix(nrow=181,ncol=0)
    mn <- matrix(nrow=181,ncol=0)
    me <- matrix(nrow=181,ncol=0)
    valNum <- 0
    for(i in (lubridate::year(as.Date(dat2$dates[1]))+1):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
      subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
      valNum <- valNum + 1
      if(valNum>length(cVals)){ ##Done to set the c and d values at the last partial year as those of the last full year
        valNum <- length(cVals)
      }
      c <- cVals[valNum]
      d <- dVals[valNum]

      p <- cbind(p,rescale(yseq=subDat$p,c=c,d=d))
      if(forecastType=="logisticCov"){
        Sf <- cbind(Sf,subDat$Sf)
      }
      mn <- cbind(mn,rescale(yseq=subDat$mn,c=c,d=d))
      me <- cbind(me,rescale(yseq=subDat$me,c=c,d=d))
    }

    ##Forecast Sf information
    if(forecastType=="logisticCov"){
      Sfmu <- matrix(nrow=nrow(mn),ncol=ncol(mn))
      Sfprec <- matrix(nrow=nrow(mn),ncol=ncol(mn))

      colNum <- ncol(Sfmu)
      rowNums <- seq((lubridate::day(as.Date(endDate))+forecastLength-length(SfsPrec)),lubridate::day(as.Date(endDate))+forecastLength,1)
      print(rowNums)
      for(i in (1:length(rowNums))){
        Sfmu[rowNums[i],colNum] <- SfsFuture[i]
        Sfprec[rowNums[i],colNum] <- SfsPrec[i]
      }
    }

    dataFinal <- list(p=p,mn=mn,me=me)
    dataFinal$n <- 181
    dataFinal$N <- ncol(dataFinal$p)
    dataFinal$x_ic <- 0
    dataFinal$tau_ic <- 1/(phenoData$g_std[1]**2)
    dataFinal$q <- lubridate::day(as.Date(endDate))+forecastLength
    if(forecastType=="logistic"){
      plot(seq(1,length(dataFinal$p[,1]),1),dataFinal$p[,1],pch=20)
      j.model <- logisticPhenoModel(data=dataFinal,nchain=nchain)
      variableNames <- c("p.proc","p.PC","p.ME","p.MN","x","r")
    }

    if(forecastType=="logisticCov"){
      dataFinal$Sf <- Sf
      dataFinal$Sfmu <- Sfmu
      dataFinal$Sfprec <- Sfprec
      dataFinal$L <- length(SfsPrec)
      print(Sfmu[28:43,7])
      print(Sfprec[28:43,7])
      j.model <- logisticCovPhenoModel(data=dataFinal,nchain=nchain)
      variableNames <- c("p.PC","p.MN","p.ME","p.proc","x","b0","b1")
    }
  }

  #   j.model <- logisticPhenoModel(data=data,nchain=nchain)
  #   variableNames <- c("p.proc","p.PC","p.ME","p.MN","x","r")
  #
  # }
  #
  # else if(forecastType=="logisticCov"){
  #   Tairs <- download_US_WCr_met(start_date=startDate,end_date=endDate)
  #   TairsForecast <- load_GEFS_Forecast("PhenologyForecastData/NOAA_GEFS.willowcreek.1.2019-01-25T12:00.2019-02-10T12:00.nc")
  #   #print(length(Tairs))
  #   #print("Done downloading met data")
  #   days <- seq(as.Date(startDate),(as.Date(endDate)+forecastLength),"day")
  #   #print(length(days))
  #   #print("Done with days")
  #   Tairs <- c(Tairs,TairsForecast)
  #   #print(length(Tairs))
  #   newMonths <- lubridate::month(days)
  #   newYears <- lubridate::year(days)
  #   #print("Done with newYears")
  #   URL <- "https://phenocam.sr.unh.edu/data/archive/willowcreek/ROI/willowcreek_DB_1000_1day.csv"
  #
  #   phenoData <- download.phenocam(URL)
  #   p.old <- phenoData$gcc_mean
  #   time.old <-  as.Date(phenoData$date)
  #   p <- rep(NA,length(days))
  #   for(i in 1:length(p.old)){
  #     p[which(days==time.old[i])] <- p.old[i]
  #   }
  #   #print("Done with p")
  #   #timeForecast <- c(days,seq.Date(from=endDate,by="day",length.out=forecastLength))
  #   dat2 <- data.frame(dates=days,years=newYears,months=newMonths,Tairs=Tairs,p=p)
  #   dat2$mn <- prepareMODIS(startDate=startDate,endDate=endDate,metric="NDVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)
  #   dat2$me <- prepareMODIS(startDate=startDate,endDate=endDate,metric="EVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)
  #   #print("Done with MODIS")
  #   #print(dim(dat2))
  #   dat2 <- dat2[dat2$months%in%seq(1,6,1),]
  #   #print(dim(dat2))
  #
  #
  #   data2 <- data.frame(p=rescaleObs(vals=dat2$p,times=dat2$dates),Sf=calSf(dat2$Tairs,days=dat2$dates))
  #   #data2$p <- rescaleObs(vals=dat2$p,times=dat2$dates)
  #   #data2$Sf <- calSf(dat2$Tairs,days=dat2$dates)
  #   data2$mn <- dat2$mn
  #   data2$me <- dat2$me
  #
  #   p <- matrix(nrow=181,ncol=0)
  #   Sf <- matrix(nrow=181,ncol=0)
  #   mn <- matrix(nrow=181,ncol=0)
  #   me <- matrix(nrow=181,ncol=0)
  #   for(i in lubridate::year(as.Date(dat2$dates[1])):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){
  #     subDat <- data2[lubridate::year(as.Date(dat2$dates))==i,]
  #     p <- cbind(p,subDat$p)
  #     Sf <- cbind(Sf,subDat$Sf)
  #     mn <- cbind(mn,subDat$mn)
  #     me <- cbind(me,subDat$me)
  #   }
  #   dataFinal <- list(p=p,Sf=Sf,mn=mn,me=me)
  #
  #
  #   dataFinal$n <- 181
  #   dataFinal$N <- ncol(dataFinal$p)
  #   dataFinal$x_ic <- 0
  #   dataFinal$tau_ic <- 1/(phenoData$g_std[1]**2)
  #   dataFinal$q <- lubridate::day(as.Date(endDate))+forecastLength
  #   print(dataFinal$q)
  #   print(dataFinal$N)
  #   print(dataFinal$n)
  #   j.model <- logisticCovPhenoModel(data=dataFinal,nchain=nchain)
  #   #variableNames <- c("p.PC","p.MN","p.ME","p.proc","x","b0","b1")
  #
  # }

  ##Run Model until convergence with large enough sample size
  out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames)
  return(out.burn)
  #return(j.model)
  #return(dataFinal)

}
