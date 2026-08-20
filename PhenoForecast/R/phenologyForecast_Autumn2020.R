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
#' @param GEFS_Files The filenames for the GEFS files
#' @param GEFS_Directory The directory where the GEFS files are located
#' @param baseTemp If CDD, the desired base temperature
#' @param index The desired index (GCC, NDVI, EVI, or all)
#' @param prevOutBurn The most recent previous forecast output
#' @param baseNum
#' @param iterSize
#' @param effSize
#' @param partialFile
#' @import rjags
#' @import runjags
#' @import coda
#' @import PhenologyBayesModeling
#' @export
phenologyForecast_Autumn <- function(forecastType,forecastLength=14,siteName,
                                     URLs,lat,long,dataDirectory,startDate,endDate,
                                     cValsPC,dValsPC,
                                     GEFS_Files="",GEFS_Directory,station="",
                                     baseTemp=NA, index="all",
                                     prevOutBurn=NA,
                                     baseNum=10000,
                                     iterSize=5000,
                                     effSize=5000, partialFile=FALSE){
  print(forecastType)
  nchain=5
  ###Download PhenoCam data and format
  phenoData <- matrix(nrow=0,ncol=32)
  for(u in 1:length(URLs)){
    print(URLs[u])
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

  months <- lubridate::month(days)
  years <- lubridate::year(days)
  p[p<0] <- 0
  p[p>1] <- 1


  pdf(paste(siteName,"_",endDate,"_DataPlots.pdf",sep=""),width=10,height=8)

  if(forecastType=="RW"){
    data$p <- rescaleObs(times=days,vals=p,partialStart=TRUE,cVals=cValsPC,dVals=dValsPC)
    data$n <- length(data$p)
    data$p[data$p<0] <- 0
    data$p[data$p>1] <- 1

    plot(days,data$p,pch=20,main="PhenoCam Data")
    abline(v=endDate,col="red")

    dev.off()

    j.model <- randomWalkPhenoModel(data=data,nchain=nchain,index = index)
    print("Done with creating the  random walk model")

    variableNames <- c("p.PC","p.proc","x")
    out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,baseNum = 5000,iterSize = 5000)

    ##Thin the data:
    out.mat <- as.matrix(out.burn$params)
    thinAmount <- round(nrow(out.mat)/5000,digits=0)
    out.burn2 <- list()
    out.burn2$params <- window(out.burn$params,thin=thinAmount)
    out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
    out.burn <- out.burn2
    return(out.burn)
  }else if(forecastType=="DownloadCovariates"){
    dat2 <- data.frame(dates=days,years=years,months=months,p=p)

    datTairs <- createTairs(lat=lat,long=long,dates=days,siteName=siteName,dataDirectory=dataDirectory,
                            endDate=(endDate+forecastLength),GEFS_Files=GEFS_Files,
                            GEFS_Directory=GEFS_Directory,forecastLength=forecastLength,station=station)
    dat2$TairMu <- datTairs$TairMu
    dat2$TairPrec<- datTairs$TairPrec

    dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(213,365),]
    nrowNum <- 365-213
    p <- matrix(nrow=nrowNum,ncol=0)
    TairMu <- matrix(nrow=nrowNum,ncol=0)

    TairPrec <- matrix(nrow=nrowNum,ncol=0)
    valNum <- 0
    days2 <- matrix(nrow=nrowNum,ncol=0)
    for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
      subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
      valNum <- valNum + 1
      if(valNum>length(cValsPC)){ ##Done to set the c and d values at the last partial year as those of the last full year (the lengths should be the same for the different cVals)
        valNum <- length(cValsPC)
      }
      c <- cValsPC[valNum]
      d <- dValsPC[valNum]
      p <- cbind(p,rescale(yseq=subDat$p,c=c,d=d))
      days2 <- cbind(days2,as.Date(subDat$dates))
      TairMu <- cbind(TairMu,subDat$TairMu)
      TairPrec <- cbind(TairPrec,subDat$TairPrec)
    }

    # dataFinal <- list(p=p)
    # dataFinal$n <- nrowNum
    # dataFinal$N <- ncol(dataFinal$p)
    # dataFinal$x_ic <- 0.99
    # dataFinal$tau_ic <- 1/(phenoData$g_std[1]**2)
    # dataFinal$q <- as.numeric(format(endDate,"%j"))+forecastLength - 213

    plot(days2,dataFinal$p,pch=20,main="PhenoCam Data")
    abline(v=endDate,col="red")


    print("Done with formating data")

    dataFinal$TairMu <- TairMu
    dataFinal$TairPrec <- TairPrec
    ##Save Tair data:
    TairData <- list(TairMu=dataFinal$TairMu,TairPrec=dataFinal$TairPrec)
    TairFileName <- paste(siteName,"_",endDate,"_meanTair.R",sep="")
    save(TairData,file=TairFileName)

    plot(days2,dataFinal$TairMu,pch=20,main="Tair")
    abline(v=endDate,col="red")
    dev.off()

    print("Done creating the CDD covariate model")


  }else if(forecastType=="general_noVariation"){


  }

  else{
    print("Forecast type not known!!!")
  }
}

