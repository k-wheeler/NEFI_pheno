#!/usr/bin/env Rscript

##Install latest versions of packages
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary") 
install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary") 

library("PhenoForecast",lib.loc="/projectnb/dietzelab/kiwheel/Rlibrary/")
library("PhenologyBayesModeling",lib.loc="/projectnb/dietzelab/kiwheel/Rlibrary/")
library("coda")
library("dplyr")
library("rjags")
library("MODISTools")
library(doParallel)


##Set and register cores for parallel
n.cores <- 10
registerDoParallel(cores=n.cores)

##Read in data
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 14

iseq <- c(1)
iseq <- c(1,2,3,4,6,15,16,20,18,24)

index <- "GCC"
#endDate <- (Sys.Date()-1)
dates <- seq(as.Date("2019-08-01"),as.Date("2019-12-10"),"week")
#print(endDate)
#Create Forecast outputs
output <- 
  foreach(i=iseq) %dopar% {
    #for(i in iseq){
    for(d in 1:length(dates)){
      endDate <- dates[d]
      print(endDate)
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
      ##Download MODIS data
      
      downloadMODIS(startDate=startDate,endDate=endDate,dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
      print("downloaded MODIS")
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
      outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_",index,"_outBurn.RData",sep="")
      if(!file.exists(outputFile)){
        outBurnRW <- phenologyForecast_Autumn(forecastType = "RW",forecastLength = forecastLength,
                                              siteName=siteName,URLs=URL,lat=lat,long=long,
                                              dataDirectory=dataDirectory,startDate,endDate,cValsPC=cMeans.p,
                                              dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,
                                              cValsME=cMeans.me,dValsME=dMeans.me,season="fall",index="GCC")
        if(typeof(outBurnRW)!=typeof(FALSE)){
          save(outBurnRW,file=outputFile)
        }
      }
      
      #Create logistic forecast if needed
      outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_basicLog_",index,"_outBurn.RData",sep="")
      if(!file.exists(outputFile)){
        outBurnL <- phenologyForecast_Autumn(forecastType = "log",forecastLength = forecastLength,
                                             siteName=siteName,URLs=URL,lat=lat,long=long,dataDirectory=dataDirectory,
                                             startDate=startDate,endDate=endDate,cValsPC=cMeans.p,dValsPC=dMeans.p,
                                             cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me,
                                             season="fall",index="GCC")
        if(typeof(outBurnL)!=typeof(FALSE)){
          save(outBurnL,file=outputFile)
        }
      }
    }
    
  }


