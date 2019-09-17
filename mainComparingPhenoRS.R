#!/usr/bin/env Rscript
#Sys.setenv(LIBCURL_BUILD="winssl")
#install.packages("https://github.com/jeroen/curl/archive/master.tar.gz", repos = NULL)

install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
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
dataDirectory <- paste(getwd(),"/",sep="")

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)

DB.vars <- c("TranF","bF","TranS","bS","c","d","prec","k")
SH.vars <- c("Tran","b","c","d","k","r","prec")

#detect cores.
#n.cores <- detectCores()
n.cores <- 5


#register the cores.
#registerDoParallel(cores=n.cores)
#i <- 14
i <- 1

iseq <- c(seq(1,6),seq(8,10),seq(15,20))
#iseq <- c(seq(9,10),seq(15,20))
#output <- 
#foreach(i=1:nrow(siteData)) %dopar% {

for(i in iseq){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  URL <- as.character(siteData$URL[i])
  PFT <- as.character(siteData$PFT[i])
  print(PFT)
  lat <- as.character(siteData$Lat[i])
  long <- as.character(siteData$Long[i])
  TZ <- as.character(siteData$TZ[i])
  downloadMODIS(startDate=as.Date("2017-07-01"),endDate=as.Date("2018-06-30"),
                dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
  
#  print(PFT)                    
  if(PFT=="DB"){
    startDay <- 182
    endDay <- 546
    xseq <- seq(startDay,endDay,1)
    fileName <- paste(siteName,"_",startDay,"_",endDay,"_PC_EIV_varBurn.RData",sep="")
    if(!file.exists(fileName)){
      j.model.PC <- createBayesModel.DB(dataSource="PC.GCC",siteName=siteName,URL=URL)#,startDay = startDay,endDay = endDay)
      PC.md.out <- runMCMC_Model(j.model=j.model.PC,variableNames = DB.vars,baseNum = 5000,iterSize = 1000,maxGBR = 5)
      save(PC.md.out,file=fileName)
    }
    ##****************
    fileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_NDVI_EIV_varBurn.RData",sep="")
    print(fileName)
    if(!file.exists(fileName)){
      j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.NDVI",siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay)
      MODIS.N.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars,baseNum=5000,iterSize = 1000)
      save(MODIS.N.md.out,file=fileName)
    }
    fileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_EVI_EIV_varBurn.RData",sep="")
    if(!file.exists(fileName)){
      j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.EVI",siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay)
      MODIS.E.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars,baseNum=5000,iterSize = 1000)
      save(MODIS.E.md.out,file=fileName)
    }
    # fileName <- paste(siteName,"_",startDay,"_",endDay,"_GOES_noon_varBurn.RData",sep="")
    # if(!file.exists(fileName)){
    # j.model.GOES <- createBayesModel.DB(dataSource="GOES.NDVI",siteName=siteName,startDay = startDay,endDay = endDay)
    # GOES.md.out <- runMCMC_Model(j.model=j.model.GOES,variableNames = DB.vars)
    # save(GOES.md.out,file=fileName)
    # }
  }
}
 
