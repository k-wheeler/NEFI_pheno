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
n.cores <- 10
registerDoParallel(cores=n.cores)

##Read in data
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 14


iseq <- c(1,2,3,4,6,15,16,20,18,24)
iseq <- c(1,2,3,6,15,18,24)
#endDate <- (Sys.Date()-1)
dates <- seq(as.Date("2019-08-01"),as.Date("2019-12-10"),"day")
#print(endDate)
#Create Forecast outputs
output <- 
  foreach(i=iseq) %dopar% {
    #for(i in iseq){
    for(d in 1:length(dates)){
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
      
      #downloadMODIS(startDate=startDate,endDate=endDate,dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
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
      
      ##Create a CDD forecast model
      #outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_NA_meanTair_outBurn_GCC.RData",sep="")
      outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_20_meanTair_outBurn_GCC.RData",sep="")
      
      if(!file.exists(outputFile)){
        partialFileName <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_20_meanTair_outBurn_GCC_PARTIAL_FILE.RData",sep="")
        outBurn <- phenologyForecast_Autumn(forecastType = "CDD_meanTair",forecastLength = forecastLength,
                                            siteName=siteName,URLs=URL,lat=lat,long=long,dataDirectory=dataDirectory,
                                            as.Date(startDate),as.Date(endDate),GEFS_Files=GEFS_files,cValsPC=cMeans.p,
                                            dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,
                                            dValsME=dMeans.me,GEFS_Directory = GEFS_Directory,station=station,
                                            season="fall",baseTemp=20, index="GCC",
                                            baseNum=20000,
                                            iterSize=10000,
                                            effSize=5000, partialFile=partialFileName)
        if(typeof(outBurn)!=typeof(FALSE)){
          save(outBurn,file=outputFile)
        }
      }
    }
  }


