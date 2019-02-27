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
#' @param cValsPC The c values for rescaling for PhenoCam
#' @param dValsPC The d values for rescaling for PhenoCam
#' @param cValsMN The c values for rescaling for MODIS NDVI
#' @param dValsMN The d values for rescaling for MODIS NDVI
#' @param cValsME The c values for rescaling for MODIS EVI
#' @param dValsME The d values for rescaling for MODIS EVI
#' @param GEFS_Files The filenames for the GEFS files
#' @param GEFS_Directory The directory where the GEFS files are located
#' @import rjags
#' @import runjags
#' @import coda
#' @import PhenologyBayesModeling
#' @export
phenologyForecast <- function(forecastType,forecastLength=16,siteName,URL,lat,long,dataDirectory,startDate,endDate,cValsPC,dValsPC,cValsMN,dValsMN,cValsME,dValsME,GEFS_Files="",GEFS_Directory){
  nchain=5
  ###Download PhenoCam data and format
  PCfileName <- paste(dataDirectory,siteName,"_",startDate,"_",endDate,"_PC_Data.RData",sep="")
  if(!file.exists(PCfileName)){
    phenoData <- download.phenocam(URL)
    save(phenoData,file=PCfileName)
  }
  load(PCfileName)
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
    data$p <- rescaleObs(times=days,vals=p,partialStart=TRUE,cVals=cValsPC,dVals=dValsPC)
    data$n <- length(data$p)
    data$mn <- rescaleObs(times=days,vals=mn,partialStart=TRUE,cVals=cValsMN,dVals=dValsMN)
    data$me <- rescaleObs(times=days,vals=me,partialStart=TRUE,cVals=cValsME,dVals=dValsME)
    print("Done with formatting data")

    j.model <- randomWalkPhenoModel(data=data,nchain=nchain)
    print("Done with creating the  random walk model")
    variableNames <- c("p.PC","p.MN","p.ME","p.proc","x")
    out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum = 5000,iterSize = 5000)
  }
  else if(forecastType=="logistic" || forecastType== "logisticCov"){
    dat2 <- data.frame(dates=days,years=years,months=months,p=p,mn=mn,me=me)
    SfsALL <- matrix(nrow=0,ncol=length(dat2$dates))
    if(forecastType=="logisticCov"){
      Tairs <- download_US_WCr_met(start_date=startDate,end_date=endDate)
      for(G in 1:length(GEFS_Files)){
        TairsForecast <- load_GEFS_Forecast(dataDirectory=GEFS_Directory,fileName=GEFS_Files[G])
        ##will need to fill in the missing (assume everything is current/up to date now)
        fileTairs <- c(Tairs,TairsForecast)
        Sfs <- calSf(Tairs=fileTairs,days=as.Date(dat2$dates))
        SfsALL <- rbind(SfsALL,Sfs)
      }
      SfsMeans <- colMeans(SfsALL)
      SfsVar <- apply(SfsALL,MARGIN=2,FUN=var)
      SfsVar[SfsVar==0] <- 0.001

      dat2$Sf <- SfsMeans
      dat2$Sfprec <- 1/SfsVar
    }

    dat2 <- dat2[dat2$months%in%seq(1,6,1),]
    p <- matrix(nrow=181,ncol=0)
    Sf <- matrix(nrow=181,ncol=0)
    mn <- matrix(nrow=181,ncol=0)
    me <- matrix(nrow=181,ncol=0)
    Sfprecs<- matrix(nrow=181,ncol=0)
    valNum <- 0
    for(i in (lubridate::year(as.Date(dat2$dates[1]))+1):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
      subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
      valNum <- valNum + 1
      if(valNum>length(cValsPC)){ ##Done to set the c and d values at the last partial year as those of the last full year (the lengths should be the same for the different cVals)
        valNum <- length(cValsPC)
      }
      c <- cValsPC[valNum]
      d <- dValsPC[valNum]

      p <- cbind(p,rescale(yseq=subDat$p,c=c,d=d))
      if(forecastType=="logisticCov"){
        Sf <- cbind(Sf,subDat$Sf)
        Sfprecs <- cbind(Sfprecs,subDat$Sfprec)
      }

      c <- cValsMN[valNum]
      d <- dValsMN[valNum]
      mn <- cbind(mn,rescale(yseq=subDat$mn,c=c,d=d))
      c <- cValsME[valNum]
      d <- dValsME[valNum]
      me <- cbind(me,rescale(yseq=subDat$me,c=c,d=d))
    }

    dataFinal <- list(p=p,mn=mn,me=me)
    dataFinal$n <- 181
    dataFinal$N <- ncol(dataFinal$p)
    dataFinal$x_ic <- 0
    dataFinal$tau_ic <- 1/(phenoData$g_std[1]**2)
    dataFinal$q <- as.numeric(format(endDate,"%j"))+forecastLength
    print("Done with formating data")
    if(forecastType=="logistic"){
      j.model <- logisticPhenoModel(data=dataFinal,nchain=nchain)
      print("Done creating the basic logistic model")
      variableNames <- c("p.proc","p.PC","p.ME","p.MN","x","r")
      out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum=10000,iterSize=5000)
    }

    if(forecastType=="logisticCov"){
      dataFinal$Sfmu <- Sf
      dataFinal$Sfprec <- Sfprecs
      j.model <- logisticCovPhenoModel(data=dataFinal,nchain=nchain)
      print("Done creating the logistic with covariate model")
      variableNames <- c("p.PC","p.MN","p.ME","p.proc","x","b1","b0")
      out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum=20000,iterSize=10000)
    }
  }

  print("Done with iterations")
  return(out.burn)

}
