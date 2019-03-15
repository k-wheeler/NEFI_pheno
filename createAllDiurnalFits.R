#!/usr/bin/env Rscript

install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("MODISTools",repo="https://cloud.r-project.org/")
#install.packages("doParallel",repo="https://cloud.r-project.org/")
library("PhenologyBayesModeling")
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

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)

s <- 9

siteName <- as.character(siteData$siteName[s])
lat <- as.numeric(siteData$Lat[s])
long <- as.numeric(siteData$Long[s])
TZ <- as.character(siteData$TZ[s])
PFT <- as.character(siteData$PFT[s])
print(siteName)    

iseq <- as.character(c(seq(1,333,1),seq(348,364,1)))

for(i in 1:length(iseq)){
  if(as.numeric(iseq[i])<10){
    iseq[i] <- paste("00",as.character(iseq[i]),sep="")
  }
  else if(as.numeric(iseq[i])<100){
    iseq[i] <- paste("0",iseq[i],sep="")
  }
}

output <- foreach(i = iseq) %dopar% {
#for(i in iseq){
#i <- iseq[4]
  if(PFT=="DB"){
  if(as.numeric(i)<182){
    year <- 2018
  }
  else{
    year <- 2017
  }
  }
  if(PFT=="SH"){
  if(as.numeric(i)<110){
    year <- 2018
  }
  else{
    year <- 2017
  }
  }
  fileName <- paste("dailyNDVI_GOES/",siteName,"_GOES_diurnal_",year,i,".csv",sep="")
  print(fileName)
  if(file.exists(fileName)){
  dat <- read.csv(fileName,header=FALSE)
  data <- list()
  print(dim(dat))
  data$x <- as.numeric(dat[3,])
  data$y <- as.numeric(dat[2,])
  if(TZ==6){
    solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(i),origin=as.Date(paste((year-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/Chicago"))$solarNoon
    solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
  }
  else if(TZ==5){
    solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(i),origin=as.Date(paste((year-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/New_York"))$solarNoon
    solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
  }
  data$mean.k <- solarNoonTime
  data$p.k <- 1
  plot(data$x,data$y,pch=20,ylim=c(0,1))
  xseq <- seq(5,20,0.01)
  #points(xseq,diurnalExp(a=0.004,c=0.2,k=12,xseq=xseq),col="red")
  
  outFileName <- paste(siteName,"_",as.character(i),"_varBurn6.RData",sep="")
  if(!file.exists(outFileName)){
    j.model <- createBayesModel.Diurnal(siteName=siteName,data)
    var.burn <- runMCMC_Model(j.model = j.model,variableNames=c("a","c","prec","k"),baseNum = 300000,iterSize=50000,maxGBR=20,maxIter=10000000)#,baseNum = 1000000,iterSize = 70000)
    if(typeof(var.burn)!=typeof(FALSE)){
      save(var.burn,file=outFileName)
    }
  }
}

}
