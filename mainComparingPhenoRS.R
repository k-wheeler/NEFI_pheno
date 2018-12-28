#!/usr/bin/env Rscript
#Sys.setenv(LIBCURL_BUILD="winssl")
#install.packages("https://github.com/jeroen/curl/archive/master.tar.gz", repos = NULL)

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("MODISTools",repo="https://cloud.r-project.org/")
#install.packages("curl",repo="https://cloud.r-project.org/")
library("PhenologyBayesModeling")
library("rjags")
library("runjags")
library("MODISTools")
library("ncdf4")
library(plyr)
library(doParallel)

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)

DB.vars <- c("TranF","bF","TranS","bS","c","d","prec","k")
SH.vars <- c("Tran","b","c","d","k","r","prec")

#detect cores.
#n.cores <- detectCores()
n.cores <- 5


#register the cores.
registerDoParallel(cores=n.cores)
#i <- 14
iseq <- c(18)
output <- 
foreach(i=1:nrow(siteData)) %dopar% {
#for(i in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  URL <- as.character(siteData$URL[i])
  PFT <- as.character(siteData$PFT[i])
  print(PFT)
  lat <- as.character(siteData$Lat[i])
  long <- as.character(siteData$Long[i])
  TZ <- as.character(siteData$TZ[i])
#  print(PFT)                    
  if(PFT=="DB"){
    startDay <- 182
    endDay <- 546
    xseq <- seq(startDay,endDay,1)
    #fileName <- paste(siteName,"_",startDay,"_",endDay,"_PC_varBurn.RData",sep="")
    #if(!file.exists(fileName)){
    #  j.model.PC <- createBayesModel.DB(dataSource="PC.GCC",siteName=siteName,URL=URL,startDay = startDay,endDay = endDay)
    #  PC.md.out <- runMCMC_Model(j.model=j.model.PC,variableNames = DB.vars,maxGBR=50)
    #  if(typeof(PC.md.out)!=typeof(FALSE)){
    #     save(PC.md.out,file=fileName)
    #  }
    #}
    ##****************
    fileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_NDVI_varBurn.RData",sep="")
    print(fileName)
    if(!file.exists(fileName)){
      j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.NDVI",siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay)
      MODIS.N.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars)
      save(MODIS.N.md.out,file=fileName)
    }
    fileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_EVI_varBurn.RData",sep="")
    if(!file.exists(fileName)){
      j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.EVI",siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay)
      MODIS.E.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars)
      save(MODIS.E.md.out,file=fileName)
    }
    #fileName <- paste(siteName,"_",startDay,"_",endDay,"_GOES_noon_varBurn.RData",sep="")
    #if(!file.exists(fileName)){
    # j.model.GOES <- createBayesModel.DB(dataSource="GOES.NDVI",siteName=siteName,startDay = startDay,endDay = endDay)
    # GOES.md.out <- runMCMC_Model(j.model=j.model.GOES,variableNames = DB.vars)
    # save(GOES.md.out,file=fileName)
    #}
  }
  else if(PFT=="SH"){
    startDay <- 110
    endDay <-90+365
    xseq <- seq(startDay,endDay,1)
    SH.vars <- c("Tran","b","c","d","k","r","prec")
    # fileName <- paste(siteName,"_PC_varBurn.RData",sep="")
    # if(!file.exists(fileName)){
    #   j.model.PC <- createBayesModel.SH(dataSource="PC.GCC",siteName=siteName,URL=URL,startDay = startDay,endDay = endDay)
    #   PC.md.out <- runMCMC_Model(j.model=j.model.PC,variableNames = SH.vars)
    #   save(PC.md.out,file=fileName)
    # }
    # fileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_NDVI_varBurn.RData",sep="")
    # if(!file.exists(fileName)){
    #   j.model.MODIS <- createBayesModel.SH(dataSource="MODIS.NDVI",siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay)
    #   MODIS.N.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = SH.vars, baseNum = 100000,iterSize = 100000)
    #   save(MODIS.N.md.out,file=fileName)
    # }
    fileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_EVI_varBurn.RData",sep="")
    if(!file.exists(fileName)){
      j.model.MODIS <- createBayesModel.SH(dataSource="MODIS.EVI",siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay)
      MODIS.E.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = SH.vars)
      save(MODIS.E.md.out,file=fileName)
    }
    #fileName <- paste(siteName,"_GOES_noon_varBurn.RData",sep="")
    #if(!file.exists(fileName)){
    # j.model.GOES <- createBayesModel.SH(dataSource="GOES.NDVI",siteName=siteName,startDay = startDay,endDay = endDay)
    # GOES.md.out <- runMCMC_Model(j.model=j.model.GOES,variableNames = SH.vars)
    # save(GOES.md.out,file=fileName)
    #}
  }
}
 
