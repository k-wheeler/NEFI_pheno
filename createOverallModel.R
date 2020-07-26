#!/usr/bin/env Rscript

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
library("ncdf4")
library(plyr)
library("PhenologyBayesModeling")
library(doParallel)
library("rjags")
library("runjags")


#detect cores.
n.cores <- 15

#register the cores.
registerDoParallel(cores=n.cores)

#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
savePath <- paste(getwd(),"/PhenologyForecastData/GOES_DiurnalFits/",sep="")
startDate=as.Date("2019-01-01")
endDate=as.Date("2019-12-31")
yr <- 2019
iseq <- c(seq(1,10),seq(11,13),15)
output <- 
  foreach(s = iseq) %dopar% {
    #for(s in iseq){
    siteName <- as.character(siteData[s,1])
    print(siteName)
    outDataFile <- paste(siteName,"_",yr,"_diurnalFitDataModel.RData",sep="")
    #if(!file.exists(outDataFile)){
    diurnalFits <- intersect(dir(path="PhenologyForecastData/GOES_DiurnalFits",pattern="varBurn"),
                             dir(path="PhenologyForecastData/GOES_DiurnalFits",
                                 pattern=paste(siteName,"_",yr,sep="")))
    c.vals <- numeric()
    prec.vals <- numeric()
    days <- numeric()
    cts <- numeric()
    for(i in 1:length(diurnalFits)){
      print(diurnalFits[i])
      load(paste("PhenologyForecastData/GOES_DiurnalFits/",diurnalFits[i],sep=""))
      if(typeof(var.burn)!=typeof(FALSE)){
        out.mat <- data.frame(as.matrix(var.burn))
        print(colnames(out.mat))
        c <- mean(out.mat$c)
        prec <- 1/var(out.mat$c)
        
        dy <- strsplit(diurnalFits[i],"_")[[1]][2]
        
        dayDataFile <- intersect(dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=siteName),
                                 dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=paste(dy,".csv",sep="")))
        #dayDataFile <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/","GOES_NDVI_Diurnal",siteName,dy,".csv",sep="")
        print(dayDataFile)
        dayData <- read.csv(paste("PhenologyForecastData/GOES_NDVI_DiurnalData/",dayDataFile,sep=""),header=FALSE)
        ct <- length(dayData[2,][!is.na(dayData[2,])])
        if(ct>1){
          c.vals <- c(c.vals,c)
          prec.vals <- c(prec.vals,prec)
          days <- c(days,substr(dy,5,7))
          cts <- c(cts,ct)
        }
      }
    }
    data <- list()
    data$x <- as.numeric(days)
    data$y <- as.numeric(c.vals)
    data$obs.prec <- as.numeric(prec.vals)
    data$n <- length(data$x)
    save(data,file=outDataFile)
    save(cts, file=paste(siteName,"_",yr,"_counts.RData",sep=""))
    print("Done with creating Data")
    varBurnFileName <- paste(siteName,"_",yr,"_GOES_NDVI_overall_varBurn.RData",sep="")
    if(!file.exists(varBurnFileName)){
      load(outDataFile)
      j.model <- createBayesModel.DB_Overall(data=data,seasonOrder="SF")
      var.burn <- runMCMC_Model(j.model=j.model,variableNames = c("TranS","bS","TranF","bF","d","c","prec"),baseNum=20000,iterSize=10000)
      save(var.burn,file=varBurnFileName)
    }
  }

