#!/usr/bin/env Rscript

##Install latest versions of packages
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary") 
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary") 

library("PhenoForecast",lib.loc="/projectnb/dietzelab/kiwheel/Rlibrary/")
library("PhenologyBayesModeling",lib.loc="/projectnb/dietzelab/kiwheel/Rlibrary/")
library("coda")
library("dplyr")
library("rjags")
library("MODISTools")
library(doParallel)


##Set and register cores for parallel
n.cores <- 15
registerDoParallel(cores=n.cores)

##Read in data
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 14

endDate <- (Sys.Date()-1)

iseq <- c(seq(1,6),8,9,10,seq(15,27))
iseq <- c(1,3,6,8,15,16,17,18,19)
#iseq <- c(1,3,6,8,18)
dates <- seq(as.Date("2019-01-23"),as.Date("2019-06-06"),"day")
dates <- seq(as.Date("2019-01-23"),as.Date("2019-02-23"),"day")

output <- foreach(d=1:length(dates)) %dopar% {
#for(d in 1:length(dates)){
  endDate <- dates[d]
  print(endDate)
  #Create Forecast outputs
  #output <- 
    #foreach(i=iseq) %dopar% {
    for(i in iseq){
      siteName <- as.character(siteData[i,1])
      print(siteName)
      GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/",siteName,"/",endDate,"/",sep="")
      GEFS_files <- dir(path=GEFS_Directory,pattern="NOAA_GEFS")
      
      URL <- as.character(siteData$URL[i])
      URL2 <- as.character(siteData$URL2[i])
      URL3 <- as.character(siteData$URL3[i])
      if(nchar(URL2)>0){
        URL <- c(URL,URL2)
        if(nchar(URL3)>0){
          URL <- c(URL,URL3)
        }
      }
      lat <- as.numeric(siteData[i,2])
      long <- as.numeric(siteData[i,3])
      startDate <- as.Date(siteData[i,7])
      station <- as.character(siteData$metStation[i])
      ##Download new MODIS data (not working in separate functions for some reason...)
      ##Download DQF file if there are no previous ones 
      #downloadMODIS(startDate=startDate,endDate=endDate,metric="rel",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
      
      #downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
      #downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
      
      ##Load rescaling data
      rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_PC.csv",sep="")
      rescaleData <- read.csv(rescaleFile,header=TRUE)
      cMeans.p <- rescaleData$cMeans.p
      dMeans.p <- rescaleData$dMeans.p
      rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_MN.csv",sep="")
      rescaleData <- read.csv(rescaleFile,header=TRUE)
      cMeans.mn <- rescaleData$cMeans.mn
      dMeans.mn <- rescaleData$dMeans.mn
      rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_ME.csv",sep="")
      rescaleData <- read.csv(rescaleFile,header=TRUE)
      cMeans.me <- rescaleData$cMeans.me
      dMeans.me <- rescaleData$dMeans.me
      
      saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
      print(saveDirectory)
      dir.create(saveDirectory)
      ##Create Random Walk forecast if needed
      outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
      if(!file.exists(outputFile)){
        outBurnRW <- phenologyForecast(forecastType = "randomWalk",forecastLength = forecastLength,siteName=siteName,URLs=URL,lat=lat,long=long,dataDirectory=dataDirectory,startDate,endDate,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me)
        if(typeof(outBurnRW)!=typeof(FALSE)){
          save(outBurnRW,file=outputFile)
        }
      }
      
      ##Create logistic forecast if needed
      outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_logistic_outBurn.RData",sep="")
      if(!file.exists(outputFile)){
        outBurnL <- phenologyForecast(forecastType = "logistic",forecastLength = forecastLength,siteName=siteName,URLs=URL,lat=lat,long=long,dataDirectory=dataDirectory,startDate=startDate,endDate=endDate,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me)
        if(typeof(outBurnL)!=typeof(FALSE)){
          save(outBurnL,file=outputFile)
        }
      }
      
      ##Create a Logistic with Covariate Model
      # outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_LC_outBurn.RData",sep="")
      # if(!file.exists(outputFile)){
      #   outBurnLC <- phenologyForecast(forecastType = "logisticCov",forecastLength = forecastLength,siteName=siteName,URLs=URL,lat=lat,long=long,dataDirectory=dataDirectory,as.Date(startDate),as.Date(endDate),GEFS_Files=GEFS_files,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me,GEFS_Directory = GEFS_Directory,station=station)
      #   if(typeof(outBurnLC)!=typeof(FALSE)){
      #     save(outBurnLC,file=outputFile)
      #   }
      # }
      outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_LC2_outBurn.RData",sep="")
      if(!file.exists(outputFile)){
       outBurnLC2 <- phenologyForecast(forecastType = "logisticCov2",forecastLength = forecastLength,siteName=siteName,URLs=URL,lat=lat,long=long,dataDirectory=dataDirectory,as.Date(startDate),as.Date(endDate),GEFS_Files=GEFS_files,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me,GEFS_Directory = GEFS_Directory,station=station)
       if(typeof(outBurnLC2)!=typeof(FALSE)){
         save(outBurnLC2,file=outputFile)
       }
      }
    }
}

