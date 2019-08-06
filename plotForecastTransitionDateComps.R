library(PhenoForecast)
library(PhenologyBayesModeling)
library(coda)
library(rjags)
library(doParallel)
library(ecoforecastR)
library("MODISTools")
library(doParallel)

##Set and register cores for parallel
n.cores <- 8
registerDoParallel(cores=n.cores)

#########Make phenofits for the last year of data
startDate <- as.Date("2018-07-01")
endDate <- as.Date("2019-06-30")
startDay <- 182
endDay <- 365+181
forecastLength <- 0
iseq <- c(1,2,3,4,6,15,16,18,20,21,22,23,24)
iseq <- c(1,2,3,4,6,15,16,20)
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
forecastDataFolder <- "PhenologyForecastData/ForecastOutputs/AllForecasts/"

#i <- 1
output <- 
  foreach(i=iseq) %dopar% {
    siteName <- as.character(siteData[i,1])
    print(siteName)
    siteStartDate <- as.character(siteData$startDate[i])
    outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_forecast_spring2019_varBurn.RData",sep="")
    
    ##First need 50% transition date estimate from the model:
    load(outFileName)
    var.mat <- data.frame(as.matrix(varBurn))
    
    
    tran50 <- quantile(var.mat$TranS,c(0.025,0.5,0.975))-365
    plotDays <- seq(23,157)
    allPlotDays <- seq(-10,170)
    allDates <- c(seq(as.Date("2019-01-23"),as.Date("2019-01-25"),"day"),
                  as.Date("2019-02-03"),as.Date("2019-02-05"),
                  seq(as.Date("2019-02-07"),as.Date("2019-06-06"),"day"))
    
    pdf(file=paste(siteName,"_PhenologyForecast_accuracy.pdf",sep=""),height=12,width=6)
    par(mfrow=c(3,1))
    #plot(plotDays,rep(NA,length(plotDays)),ylim=c(1,160),xlim=c(1,160),xlab="Forecast Day",ylab="Transition Estimate")
    ###Random Walk
    plot(plotDays,rep(NA,length(plotDays)),ylim=c(100,160),xlim=c(1,160),
         xlab="Forecast Day",ylab="Transition Estimate", main="Random Walk")
    
    ciEnvelope(allPlotDays,rep(tran50[3],length(allPlotDays)),rep(tran50[1],length(allPlotDays)),col="gray")
    abline(h=tran50[2],col="black",lwd=2)
    abline(v=tran50[2],col="black",lwd=2)
    
    allDates2 <- as.Date("2019-01-23")
    minsCI <- matrix(ncol=0,nrow=3)
    for(d in 1:length(allDates)){
      #for(d in 1:10){
      calEndDate <- allDates[d]
      #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
      #calEndDate <- as.Date("2019-05-20")
      print(calEndDate)
      forStartDate <- calEndDate + 1
      forEndDate <- (forStartDate+13)
      forDates <- seq(forStartDate,forEndDate,"day")
      yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
      
      forecastLength <- 14
      plotDates <- as.numeric(format(forDates,"%j"))
      RWfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_randomWalk_outBurn.RData",sep="")
      if(file.exists(RWfileName)){
        #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
        load(RWfileName)
        out.mat.par <- data.frame(as.matrix(outBurnRW$params))
        dayNumber <- dim(as.matrix(outBurnRW$predict))[2]-as.numeric(format(forEndDate,"%j"))+1
        out.mat <- as.matrix(outBurnRW$predict)
        out.mat <- out.mat[,dayNumber:ncol(out.mat)]
        mins <- numeric()
        for(s in 1:nrow(out.mat)){
          newVals <- which(out.mat[s,]>0.50)
          newMin <- min(newVals)
          mins <- c(mins,newMin)
        }
        minsCI <- cbind(minsCI,quantile(mins,c(0.025,0.50,0.975),na.rm=TRUE))
        allDates2 <- c(allDates2,calEndDate)
      }
    }
    allDates3 <- allDates2[2:length(allDates2)]
    ciEnvelope(format(allDates3,"%j"),minsCI[1,],minsCI[3,],col="lightblue")
    lines(format(allDates3,"%j"),minsCI[2,],col="darkblue",lwd=2)
    
    ###Logistic
    plot(plotDays,rep(NA,length(plotDays)),ylim=c(100,160),xlim=c(1,160),
         xlab="Forecast Day",ylab="Transition Estimate", main="Logistic")
    
    ciEnvelope(allPlotDays,rep(tran50[3],length(allPlotDays)),rep(tran50[1],length(allPlotDays)),col="gray")
    abline(h=tran50[2],col="black",lwd=2)
    abline(v=tran50[2],col="black",lwd=2)
    
    allDates2 <- as.Date("2019-01-23")
    minsCI <- matrix(ncol=0,nrow=3)
    for(d in 1:length(allDates)){
      #for(d in 1:10){
      calEndDate <- allDates[d]
      #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
      #calEndDate <- as.Date("2019-05-20")
      print(calEndDate)
      forStartDate <- calEndDate + 1
      forEndDate <- (forStartDate+13)
      forDates <- seq(forStartDate,forEndDate,"day")
      yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
      
      forecastLength <- 14
      plotDates <- as.numeric(format(forDates,"%j"))
      LfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_logistic_outBurn.RData",sep="")
      if(file.exists(LfileName)){
        #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
        load(LfileName)
        out.mat.par <- data.frame(as.matrix(outBurnL$params))
        dayNumber <- dim(as.matrix(outBurnL$predict))[2]-as.numeric(format(forEndDate,"%j"))+1
        out.mat <- as.matrix(outBurnL$predict)
        out.mat <- out.mat[,dayNumber:ncol(out.mat)]
        mins <- numeric()
        for(s in 1:nrow(out.mat)){
          newVals <- which(out.mat[s,]>0.50)
          newMin <- min(newVals)
          mins <- c(mins,newMin)
        }
        minsCI <- cbind(minsCI,quantile(mins,c(0.025,0.50,0.975),na.rm=TRUE))
        allDates2 <- c(allDates2,calEndDate)
      }
    }
    allDates3 <- allDates2[2:length(allDates2)]
    ciEnvelope(format(allDates3,"%j"),minsCI[1,],minsCI[3,],col="lightblue")
    lines(format(allDates3,"%j"),minsCI[2,],col="darkblue",lwd=2)
    
    
    ###Logistic Covariate
    plot(plotDays,rep(NA,length(plotDays)),ylim=c(100,160),xlim=c(1,160),xlab="Forecast Day",ylab="Transition Estimate",main="Logistic Covariate")
    ciEnvelope(allPlotDays,rep(tran50[3],length(allPlotDays)),rep(tran50[1],length(allPlotDays)),col="gray")
    abline(h=tran50[2],col="black",lwd=2)
    abline(v=tran50[2],col="black",lwd=2)
    
    allDates2 <- as.Date("2019-01-23")
    minsCI <- matrix(ncol=0,nrow=3)
    for(d in 1:length(allDates)){
      #for(d in 1:10){
      calEndDate <- allDates[d]
      #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
      #calEndDate <- as.Date("2019-05-20")
      print(calEndDate)
      forStartDate <- calEndDate + 1
      forEndDate <- (forStartDate+13)
      forDates <- seq(forStartDate,forEndDate,"day")
      yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
      
      forecastLength <- 14
      plotDates <- as.numeric(format(forDates,"%j"))
      LCfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
      if(file.exists(LCfileName)){
        #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
        load(LCfileName)
        out.mat.par <- data.frame(as.matrix(outBurnLC2$params))
        dayNumber <- dim(as.matrix(outBurnLC2$predict))[2]-as.numeric(format(forEndDate,"%j"))+1
        out.mat <- as.matrix(outBurnLC2$predict)
        out.mat <- out.mat[,dayNumber:ncol(out.mat)]
        mins <- numeric()
        for(s in 1:nrow(out.mat)){
          newVals <- which(out.mat[s,]>0.50)
          newMin <- min(newVals)
          mins <- c(mins,newMin)
        }
        minsCI <- cbind(minsCI,quantile(mins,c(0.025,0.50,0.975),na.rm=TRUE))
        allDates2 <- c(allDates2,calEndDate)
      }
    }
    allDates3 <- allDates2[2:length(allDates2)]
    ciEnvelope(format(allDates3,"%j"),minsCI[1,],minsCI[3,],col="lightblue")
    lines(format(allDates3,"%j"),minsCI[2,],col="darkblue",lwd=2)
    
    dev.off()
  }