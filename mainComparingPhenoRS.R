#!/usr/bin/env Rscript

#install.packages("https://github.com/jeroen/curl/archive/master.tar.gz", repos = NULL)

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("MODISTools",repo="https://cloud.r-project.org/")
#install.packages("curl",repo="https://cloud.r-project.org/")
library("PhenologyBayesModeling")
library("PhenoForecast")
library("rjags")
library("runjags")
library("MODISTools")
library("ncdf4")
library(plyr)
library(doParallel)

#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)


#detect cores.
#n.cores <- detectCores()
n.cores <- 16

#register the cores.
registerDoParallel(cores=n.cores)

createAllPhenoFits <- function(startDate,endDate,siteData){
  DB.vars <- c("TranF","bF","TranS","bS","c","d","prec")
  dataDirectory <- paste(getwd(),"/",sep="")
  iseq <- 1:nrow(siteData)
  output <- 
  foreach(i=1:nrow(siteData)) %dopar% {
  #for(i in iseq){
    siteName <- as.character(siteData$siteName[i])
    print(siteName)
    URL <- as.character(siteData$URL[i])
    PFT <- as.character(siteData$PFT[i])
    print(PFT)
    lat <- as.character(siteData$Lat[i])
    long <- as.character(siteData$Long[i])
    TZ <- as.character(siteData$TZ[i])
    
    fileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_EIV_varBurn_SF.RData",sep="")
    if(!file.exists(fileName)){
      j.model.PC <- createBayesModel.DB(dataSource="PC.GCC",siteName=siteName,URL=URL,seasonOrder = "SF",startDate = startDate,
                                        endDate=endDate)
      PC.md.out <- runMCMC_Model(j.model=j.model.PC,variableNames = DB.vars,baseNum = 5000,iterSize = 5000,maxGBR = 5)
      save(PC.md.out,file=fileName)
    }
    ##****************
    metric <- "NDVI"
    fileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_NDVI_EIV_varBurn.RData",sep="")
    print(fileName)
    if(!file.exists(fileName)){
      j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.NDVI",siteName=siteName,lat=lat,long=long,
                                           startDate = startDate,endDate = endDate,seasonOrder = "SF")
      MODIS.N.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars,baseNum=5000,iterSize = 5000)
      save(MODIS.N.md.out,file=fileName)
    }
    fileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_EVI_EIV_varBurn.RData",sep="")
    print(fileName)
    metric <- "EVI"
    if(!file.exists(fileName)){
      print("MODIS File Downloaded")
      j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.EVI",siteName=siteName,lat=lat,long=long,
                                           startDate = startDate,endDate = endDate,seasonOrder = "SF")
      MODIS.E.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars,baseNum=5000,iterSize = 5000)
      save(MODIS.E.md.out,file=fileName)
    }
  }
  
}
createAllPhenoFits(startDate=as.Date("2018-01-01"),endDate=as.Date("2018-12-31"),siteData=siteData)
createAllPhenoFits(startDate=as.Date("2019-01-01"),endDate=as.Date("2019-12-31"),siteData=siteData)
