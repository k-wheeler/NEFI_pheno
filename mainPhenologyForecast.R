library("PhenoForecast")
library("coda")
library(doParallel)

##Set and register cores for parallel
n.cores <- 6
registerDoParallel(cores=n.cores)

##Read in data
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
directory="PhenologyForecastData/"
forecastLength <- 100

endDate <- (Sys.Date()-1)

#Create Forecast outputs
output <- 
  foreach(i=1:nrow(siteData)) %dopar% {
    siteName <- as.character(siteData[i,1])
    print(siteName)
    URL <- as.character(siteData[i,4])
    lat <- as.numeric(siteData[i,2])
    long <- as.numeric(siteData[i,3])
    startDate <- as.Date(siteData[i,7])
    ##Download new MODIS data
    downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=directory,lat=lat,long=long,siteName=siteName)
    downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=directory,lat=lat,long=long,siteName=siteName)
    
    ##Create Random Walk forecast if needed
    outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
    if(!file.exists(outputFile)){
      outBurn <- phenologyForecast(forecastType = "randomWalk",forecastLength = forecastLength,siteName=siteName,URL=URL,lat=lat,long=long,dataDirectory=directory,startDate,endDate)
      if(typeof(outBurn)!=typeof(FALSE)){
        save(outBurn,file=outputFile)
      }
    }
    
    ##Create logistic forecast if needed
    outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_logistic_outBurn.RData",sep="")
    if(!file.exists(outputFile)){
      outBurn <- phenologyForecast(forecastType = "logistic",forecastLength = forecastLength,siteName=siteName,URL=URL,lat=lat,long=long,dataDirectory=directory,startDate,endDate)
      if(typeof(outBurn)!=typeof(FALSE)){
        save(outBurn,file=outputFile)
      }
    }
  }

##Plotting Forecasts
# pdf(file=paste("ForecastOutputs_",endDate,".pdf",sep=""),width=12,height=8)
# for(i in 1:16){
#   siteName <- as.character(siteData[i,1])
#   print(siteName)
#   URL <- as.character(siteData[i,4])
#   lat <- as.numeric(siteData[i,2])
#   long <- as.numeric(siteData[i,3])
#   startDate <- as.Date(siteData[i,7])
#   outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
#   load(outputFile)
#   plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,dataSource="PhenoCam",outBurn=outBurn,forecastType = "randomWalk")
#   
#   outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_logistic_outBurn.RData",sep="")
#   load(outputFile)
#   plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,dataSource="PhenoCam",outBurn=outBurn,forecastType = "logistic")
# }
# dev.off()
