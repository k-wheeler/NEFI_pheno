#!/usr/bin/env Rscript
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",
#                 repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary")
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/GOESDiurnalNDVI",
#                 repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary")

library("PhenologyBayesModeling")
library("GOESDiurnalNDVI")
library("rjags")
library("runjags")
#library("MODISTools")
library("doParallel")
library("suncalc")

#detect cores.
#n.cores <- detectCores()
n.cores <- 10

#register the cores.
registerDoParallel(cores=n.cores)

#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
savePath <- paste(getwd(),"/PhenologyForecastData/GOES_DiurnalFits/",sep="")

##Site Data: 
s <- 1
siteName <- as.character(siteData$siteName[s])
lat <- as.numeric(siteData$Lat[s])
long <- as.numeric(siteData$Long[s])
TZ <- as.character(siteData$TZ[s])
PFT <- as.character(siteData$PFT[s])
print(siteName)    

days <- seq(1,365)
year <- 2019

for(i in 1:length(days)){
  if(as.numeric(days[i])<10){
    days[i] <- paste("00",as.character(days[i]),sep="")
  }
  else if(as.numeric(days[i])<100){
    days[i] <- paste("0",days[i],sep="")
  }
}

#output <- foreach(day = days) %dopar% {
for(day in days){
  #fileName <- paste("dailyNDVI_GOES/",siteName,"_GOES_diurnal_",year,i,".csv",sep="")
  fileName <- paste("PhenologyForecastData/GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="")
  print(fileName)
  if(file.exists(fileName)){
    GOESdat <- read.csv(fileName,header=FALSE)
    plot(as.numeric(GOESdat[3,]),as.numeric(GOESdat[2,]),pch=20,xlim=c(0,20),ylim=c(0,1),
         ylab="NDVI",xlab="Time (Hour)")
    if(ncol(GOESdat)>1){ #Need to check that it is not an empty file
      data <- list(x=as.numeric(GOESdat[3,]),y=as.numeric(GOESdat[2,]))
      modelFitFileName <- paste(savePath,siteName,"_",year,day,"_varBurn.RData",sep="")
      if(!file.exists(modelFitFileName)){
        j.model <- createDiurnalModel(siteName=siteName,data=data)
        var.burn <- runMCMC_Model(j.model=j.model,variableNames=c("a","c","k","prec"),
                                  baseNum=40000,iterSize =30000,maxGBR = 1.80) #The baseNum and iterSize can be increased/decreased to make the code run faster if you know it will converge easier
        ##Thin:
        if(typeof(var.burn)!=typeof(FALSE)){
          out.mat <- as.matrix(var.burn)
          thinAmount <- round(nrow(out.mat)/5000,digits=0)
          var.burn <- window(var.burn,thin=thinAmount)
          save(var.burn,file=modelFitFileName)
        }
      }
    }
  }
}
  
  # dat <- read.csv(fileName,header=FALSE)
  # data <- list()
  # print(dim(dat))
  # data$x <- as.numeric(dat[3,])
  # data$y <- as.numeric(dat[2,])
  # if(TZ==6){
  #   solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(i),origin=as.Date(paste((year-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/Chicago"))$solarNoon
  #   solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
  # }
  # else if(TZ==5){
  #   solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(i),origin=as.Date(paste((year-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/New_York"))$solarNoon
  #   solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
  # }
  # data$mean.k <- solarNoonTime
  # data$p.k <- 1
  # plot(data$x,data$y,pch=20,ylim=c(0,1))
  # xseq <- seq(5,20,0.01)
  # #points(xseq,diurnalExp(a=0.004,c=0.2,k=12,xseq=xseq),col="red")
  # 
  # outFileName <- paste(siteName,"_",as.character(i),"_varBurn6.RData",sep="")
  # if(!file.exists(outFileName)){
  #   j.model <- createBayesModel.Diurnal(siteName=siteName,data)
    # var.burn <- runMCMC_Model(j.model = j.model,variableNames=c("a","c","prec","k"),baseNum = 300000,iterSize=50000,maxGBR=20,maxIter=10000000)#,baseNum = 1000000,iterSize = 70000)
    # if(typeof(var.burn)!=typeof(FALSE)){
    #   save(var.burn,file=outFileName)
    # }
  # }

