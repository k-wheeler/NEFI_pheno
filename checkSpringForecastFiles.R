#!/usr/bin/env Rscript

###Script to check for forecasts that need to be rerun
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"
forecastLength <- 14
#siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

##Forecast dates starts 2019-01-23 through 2019-06-06
dates <- seq(as.Date("2019-01-23"),as.Date("2019-06-06"),"day")
checkOutputFiles <- function(outStr){
  missingCt <- 0
  totalCt <- 0
  for(s in 1:nrow(siteData)){
    PFT <- as.character(siteData$PFT[s])
    if(PFT=="DB"){
      siteName <- as.character(siteData$siteName[s])
      #print(siteName)
      for(d in 1:length(dates)){
        outputDir <- paste(dataDirectory,"/ForecastOutputs/",siteName,"/",dates[d],"/",sep="")
        ##Check Random Walk:
        #outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
        if(length(dir(path=outputDir,pattern=outStr))==0){
          missingCt <- missingCt + 1
        }
        totalCt <- totalCt + 1
        
      }
    }
  }
  print(outStr)
  print(missingCt)
  print(missingCt/totalCt)
}
checkOutputFiles("_randomWalk_outBurn.RData")
checkOutputFiles("_logistic_outBurn.RData")
checkOutputFiles("_LC2_outBurn.RData")


