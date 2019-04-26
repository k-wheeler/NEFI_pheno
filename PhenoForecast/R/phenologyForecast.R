#' General phenology forecast function
#'
#' @param forecastType The type of forecast (randomWalk or logistic)
#' @param forecastLength The length of the forecast into the future in days
#' @param siteName The site name
#' @param URLs The URL where the site's phenocam data is located
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
phenologyForecast <- function(forecastType,forecastLength=14,siteName,URLs,lat,long,dataDirectory,startDate,endDate,cValsPC,dValsPC,cValsMN,dValsMN,cValsME,dValsME,GEFS_Files="",GEFS_Directory,station=""){
  print(forecastType)
  nchain=5
  ###Download PhenoCam data and format
  #PCfileName <- paste(dataDirectory,siteName,"_",startDate,"_",endDate,"_PC_Data.RData",sep="")
  #if(!file.exists(PCfileName)){
  phenoData <- matrix(nrow=0,ncol=32)
  for(u in 1:length(URLs)){
    phenoDataSub <- download.phenocam(URLs[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }
  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  phenoData <- phenoData[phenoData$date<endDate,]
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

  pdf(paste(siteName,"_",endDate,"_DataPlots.pdf",sep=""),width=10,height=8)

  if(forecastType=="randomWalk"){
    data$p <- rescaleObs(times=days,vals=p,partialStart=TRUE,cVals=cValsPC,dVals=dValsPC)
    data$n <- length(data$p)
    data$mn <- rescaleObs(times=days,vals=mn,partialStart=TRUE,cVals=cValsMN,dVals=dValsMN)
    data$me <- rescaleObs(times=days,vals=me,partialStart=TRUE,cVals=cValsME,dVals=dValsME)
    plot(days,data$p,pch=20,main="PhenoCam Data")
    abline(v=endDate,col="red")
    plot(days,data$mn,pch=20,main="MODIS NDVI Data")
    abline(v=endDate,col="red")
    plot(days,data$me,pch=20,main="MODIS EVI Data")
    abline(v=endDate,col="red")
    print("Done with formatting data")
    dev.off()

    j.model <- randomWalkPhenoModel(data=data,nchain=nchain)
    print("Done with creating the  random walk model")
    variableNames <- c("p.PC","p.MN","p.ME","p.proc","x")
    #out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum = 5000,iterSize = 5000)
  }else if(forecastType=="logistic" || forecastType== "logisticCov" || forecastType== "logisticCov2"|| forecastType== "logisticCov3"){
    dat2 <- data.frame(dates=days,years=years,months=months,p=p,mn=mn,me=me)
    if(forecastType=="logisticCov" || forecastType== "logisticCov2"|| forecastType== "logisticCov3"){
      datSf <- createSf(lat=lat,long=long,dates=days,siteName=siteName,dataDirectory=dataDirectory,endDate=(endDate+forecastLength),GEFS_Files=GEFS_Files,GEFS_Directory=GEFS_Directory,forecastLength=forecastLength,station=station)

      dat2$Sf <- datSf$Sf
      dat2$Sfprec <- datSf$Sfprec
    }

    #dat2 <- dat2[dat2$months%in%seq(1,6,1),]
    dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(1,181),]
    p <- matrix(nrow=181,ncol=0)
    Sf <- matrix(nrow=181,ncol=0)
    mn <- matrix(nrow=181,ncol=0)
    me <- matrix(nrow=181,ncol=0)
    Sfprecs<- matrix(nrow=181,ncol=0)
    valNum <- 0
    days2 <- matrix(nrow=181,ncol=0)
    for(i in (lubridate::year(as.Date(dat2$dates[1]))+1):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
      subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
      valNum <- valNum + 1
      if(valNum>length(cValsPC)){ ##Done to set the c and d values at the last partial year as those of the last full year (the lengths should be the same for the different cVals)
        valNum <- length(cValsPC)
      }
      c <- cValsPC[valNum]
      d <- dValsPC[valNum]

      p <- cbind(p,rescale(yseq=subDat$p,c=c,d=d))
      #print(subDat$dates)
      days2 <- cbind(days2,as.Date(subDat$dates))
      if(forecastType=="logisticCov" || forecastType== "logisticCov2" || forecastType== "logisticCov3"){
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

    plot(days2,dataFinal$p,pch=20,main="PhenoCam Data")
    abline(v=endDate,col="red")
    plot(days2,dataFinal$mn,pch=20,main="MODIS NDVI Data")
    abline(v=endDate,col="red")
    plot(days2,dataFinal$me,pch=20,main="MODIS EVI Data")
    abline(v=endDate,col="red")

    print("Done with formating data")
    if(forecastType=="logistic"){
      dev.off()
      j.model <- logisticPhenoModel(data=dataFinal,nchain=nchain)
      print("Done creating the basic logistic model")
      variableNames <- c("p.proc","p.PC","p.ME","p.MN","x","r")
      #out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum=10000,iterSize=5000)
    }else if(forecastType=="logisticCov"){
      dataFinal$Sfmu <- Sf
      dataFinal$Sfprec <- Sfprecs
      plot(days2,dataFinal$Sfmu,pch=20,main="Sf")
      abline(v=endDate,col="red")
      dev.off()
      j.model <- logisticCovPhenoModel(data=dataFinal,nchain=nchain)
      print("Done creating the logistic with covariate model")
      variableNames <- c("p.PC","p.MN","p.ME","x","p.proc","b1","b0")
      print(variableNames)
      #out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum=20000,iterSize=10000)
    }else if(forecastType=="logisticCov2"){
      dataFinal$Sfmu <- Sf
      dataFinal$Sfprec <- Sfprecs
      plot(days2,dataFinal$Sfmu,pch=20,main="Sf")
      abline(v=endDate,col="red")
      dev.off()
      print("Creating logistic with covariate model 2")
      j.model <- logisticCovPhenoModel2(data=dataFinal,nchain=nchain)
      print("Done creating the logistic with covariate model 2")
      variableNames <- c("p.PC","p.MN","p.ME","x","p.proc","b1","trans")
      print(variableNames)
      #out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum=5000,iterSize=2000)
    }else if(forecastType=="logisticCov3"){
      dataFinal$Sfmu <- Sf
      dataFinal$Sfprec <- Sfprecs
      plot(days2,dataFinal$Sfmu,pch=20,main="Sf")
      abline(v=endDate,col="red")
      dev.off()
      print("Creating logistic with covariate model 3")
      j.model <- logisticCovPhenoModel3(data=dataFinal,nchain=nchain)
      print("Done creating the logistic with covariate model 3")
      variableNames <- c("p.PC","p.MN","p.ME","x","p.proc","b1","trans")
      print(variableNames)
      #out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum=5000,iterSize=2000)
    }else{
      print("Forecast type not known!!!")
    }

  }
  print("Done with iterations")

  ##Thin the data:
  #out.mat <- as.matrix(out.burn)
  #out.burn <- window(out.burn,thin=(nrow(out.mat)/5000))

  #print("Done thinning")
  return(out.burn)
  return(FALSE)

}
