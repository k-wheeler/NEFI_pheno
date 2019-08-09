##Script to make the comprehensive plot/table comparing the predictability of the different models for different sites and different transitions
##Need to start with one site and one model and calculate how far off the transition date estimate occurs within the actual transition date CI

library(PhenoForecast)
library(PhenologyBayesModeling)
library(coda)
library(rjags)
library(doParallel)
library(ecoforecastR)
library("MODISTools")
library(doParallel)

##Set and register cores for parallel
n.cores <- 2
registerDoParallel(cores=n.cores)

#########Make phenofits for the last year of data
startDate <- as.Date("2018-07-01")
endDate <- as.Date("2019-06-30")
startDay <- 182
endDay <- 365+181
forecastLength <- 0
iseq <- c(1,2,3,4,6,15,16,18,20,21,22,23,24)
iseq <- c(1,2,3,4,6,15,16,20)
iseq <- c(18,24)
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
forecastDataFolder <- "PhenologyForecastData/ForecastOutputs/AllForecasts/"
dataFolder <- "PhenologyForecastData/"
outTable <- matrix(nrow=length(iseq),ncol=3)
i <- 1
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
    ##Random Walk
    outFileName <- paste(dataFolder, siteName, "_RW_spring50TranComp.RData",sep="")
    if(!file.exists(outFileName)){
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
      RWdat <- list(dts=allDates3,minsCI=minsCI)
      save(RWdat,file = outFileName)
    }
    
    ###Logistic
    
    outFileName <- paste(dataFolder, siteName, "_L_spring50TranComp.RData",sep="")
    if(!file.exists(outFileName)){
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
      Ldat <- list(dts=allDates3,minsCI=minsCI)
      save(Ldat,file = outFileName)
    }
    
    
    ###Logistic Covariate
    outFileName <- paste(dataFolder, siteName, "_LC_spring50TranComp.RData",sep="")
    if(!file.exists(outFileName)){
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
      LCdat <- list(dts=allDates3,minsCI=minsCI)
      save(LCdat,file = outFileName)
    
    }
    # load(outFileName)
    # potInd <- which(minsCI[1,]<tran50[1])
    # 
    # noFrstPot <- TRUE
    # j <- 1
    # while(noFrstPot){ ##Only counts it as the first day if the next three days are also included 
    #   print(j)
    #   if(potInd[j+1]==potInd[j]+1 &&
    #      potInd[j+2]==potInd[j]+2 &&
    #      potInd[j+3]==potInd[j]+3){
    #     frstDy <- potInd[j]
    #     noFrstPot <- FALSE
    #   }else{
    #     j <-  j + 1
    #   }
    #   if(j>(length(potInd)-3)){
    #     noFrstPot <- FALSE
    #     frstDy <- NA
    #   }
    # }
  }
}
