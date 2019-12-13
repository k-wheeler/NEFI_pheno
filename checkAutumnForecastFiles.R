#!/usr/bin/env Rscript

##Check for Random Walk GCC files
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"
iseq <- c(1,2,3,4,6,15,16,20,18,24)

index <- "GCC"
dates <- seq(as.Date("2019-08-01"),as.Date("2019-12-10"),"day")
print("randomWalk_GCC")
totalFiles <- 0
currentFiles <- 0
for(i in iseq){
  totalFilesSite <- 0
  currentFilesSite <- 0
  siteName <- as.character(siteData[i,1])
  startDate <- as.Date(siteData[i,7])
  print(siteName)
  for(d in 1:length(dates)){
    endDate <- dates[d]
    saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
    outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_",index,"_outBurn.RData",sep="")
    if(file.exists(outputFile)){
      currentFiles <- currentFiles + 1
      currentFilesSite <- currentFilesSite + 1
    }
    totalFiles <- totalFiles + 1
    totalFilesSite <- totalFilesSite + 1
  }
  print(currentFilesSite/totalFilesSite)
}
print("Total")
print(currentFiles/totalFiles)

print("basicLog_GCC")
totalFiles <- 0
currentFiles <- 0
for(i in iseq){
  totalFilesSite <- 0
  currentFilesSite <- 0
  siteName <- as.character(siteData[i,1])
  startDate <- as.Date(siteData[i,7])
  print(siteName)
  for(d in 1:length(dates)){
    endDate <- dates[d]
    saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
    outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_basicLog_",index,"_outBurn.RData",sep="")
    if(file.exists(outputFile)){
      currentFiles <- currentFiles + 1
      currentFilesSite <- currentFilesSite + 1
    }
    totalFiles <- totalFiles + 1
    totalFilesSite <- totalFilesSite + 1
  }
  print(currentFilesSite/totalFilesSite)
}
print("Total")
print(currentFiles/totalFiles)

print("CDD20_meanTair_GCC")
totalFiles <- 0
currentFiles <- 0
for(i in iseq){
  totalFilesSite <- 0
  currentFilesSite <- 0
  siteName <- as.character(siteData[i,1])
  startDate <- as.Date(siteData[i,7])
  print(siteName)
  for(d in 1:length(dates)){
    endDate <- dates[d]
    saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
    outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_20_meanTair_outBurn_GCC.RData",sep="")
    if(file.exists(outputFile)){
      currentFiles <- currentFiles + 1
      currentFilesSite <- currentFilesSite + 1
    }
    totalFiles <- totalFiles + 1
    totalFilesSite <- totalFilesSite + 1
  }
  print(currentFilesSite/totalFilesSite)
}
print("Total")
print(currentFiles/totalFiles)

print("CDD10_meanTair_GCC")
totalFiles <- 0
currentFiles <- 0
for(i in iseq){
  totalFilesSite <- 0
  currentFilesSite <- 0
  siteName <- as.character(siteData[i,1])
  startDate <- as.Date(siteData[i,7])
  print(siteName)
  for(d in 1:length(dates)){
    endDate <- dates[d]
    saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
    outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_10_meanTair_outBurn_GCC.RData",sep="")
    if(file.exists(outputFile)){
      currentFiles <- currentFiles + 1
      currentFilesSite <- currentFilesSite + 1
    }
    totalFiles <- totalFiles + 1
    totalFilesSite <- totalFilesSite + 1
  }
  print(currentFilesSite/totalFilesSite)
}
print("Total")
print(currentFiles/totalFiles)

print("CDD05_meanTair_GCC")
totalFiles <- 0
currentFiles <- 0
for(i in iseq){
  totalFilesSite <- 0
  currentFilesSite <- 0
  siteName <- as.character(siteData[i,1])
  startDate <- as.Date(siteData[i,7])
  print(siteName)
  for(d in 1:length(dates)){
    endDate <- dates[d]
    saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
    outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_05_meanTair_outBurn_GCC.RData",sep="")
    if(file.exists(outputFile)){
      currentFiles <- currentFiles + 1
      currentFilesSite <- currentFilesSite + 1
    }
    totalFiles <- totalFiles + 1
    totalFilesSite <- totalFilesSite + 1
  }
  print(currentFilesSite/totalFilesSite)
}
print("Total")
print(currentFiles/totalFiles)

print("CDD00_meanTair_GCC")
totalFiles <- 0
currentFiles <- 0
for(i in iseq){
  totalFilesSite <- 0
  currentFilesSite <- 0
  siteName <- as.character(siteData[i,1])
  startDate <- as.Date(siteData[i,7])
  print(siteName)
  for(d in 1:length(dates)){
    endDate <- dates[d]
    saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
    outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_00_meanTair_outBurn_GCC.RData",sep="")
    if(file.exists(outputFile)){
      currentFiles <- currentFiles + 1
      currentFilesSite <- currentFilesSite + 1
    }
    totalFiles <- totalFiles + 1
    totalFilesSite <- totalFilesSite + 1
  }
  print(currentFilesSite/totalFilesSite)
}
print("Total")
print(currentFiles/totalFiles)

  

  