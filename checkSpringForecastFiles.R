#!/usr/bin/env Rscript

###Script to check for forecasts that need to be rerun
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"
forecastLength <- 14
iseq <- c(3,6,18,15,16,4,2)
iseq <- c(1)
iseq <- c(20,21,22,23,24)
iseq <- c(2,3,4,6)#,15,16,18)#,20,21,22,23,24)
iseq <- c(20,21,22,23,24)
#iseq <- c(20)
#iseq <- c(20,21)
#iseq <- c(1,3,6,8,9,15,16,17,18,19)
#siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

##Forecast dates starts 2019-01-23 through 2019-06-06
#dates <- seq(as.Date("2019-01-23"),as.Date("2019-06-06"),"day")
dates <- c(seq(as.Date("2019-01-23"),as.Date("2019-01-25"),"day"),
           as.Date("2019-02-03"),as.Date("2019-02-05"),
	  seq(as.Date("2019-02-07"),as.Date("2019-06-06"),"day"))
        
checkOutputFiles <- function(outStr){
  missingCt <- 0
  totalCt <- 0
  for(s in iseq){
    PFT <- as.character(siteData$PFT[s])
    if(PFT=="DB"){
      siteName <- as.character(siteData$siteName[s])
      #print(siteName)
      for(d in 1:length(dates)){
        outputDir <- paste(dataDirectory,"/ForecastOutputs/",siteName,"/",dates[d],"/",sep="")
        ##Check Random Walk:
        #outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
        if(length(dir(path=outputDir,pattern=outStr))==0){
          print(paste(siteName,dates[d]))
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
#checkOutputFiles("_LC2_outBurn.RData")


