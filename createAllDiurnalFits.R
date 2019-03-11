#!/usr/bin/env Rscript

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
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
n.cores <- 6

#register the cores.
registerDoParallel(cores=n.cores)

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)

s <- 9

siteName <- as.character(siteData$siteName[s])
lat <- as.numeric(siteData$Lat[s])
long <- as.numeric(siteData$Long[s])
TZ <- as.character(siteData$TZ[s])
print(siteName)    
#diurnal.files <- dir(path="dailyNDVI_GOES",pattern=paste("GOES_Diurnal_",siteName,sep=""))
#iseq <- c(186,191,198,230,248,250,252,285)
#iseq <- c(seq(186,193),seq(195,201),206,207,211,217,230,231,seq(233,236),seq(244,254),258,259,seq(277,287),seq(297,299),seq(301,304),seq(313,315))
#iseq <- c(seq(182,321,1),seq(348,364,1))
#russellSage: #iseq <- c(seq(182,193),seq(195,203),seq(206,207),seq(211,213),seq(215,217),224,seq(227,231),seq(233,236),seq(244,250),seq(251,254),seq(258,260),seq(262,268),seq(271,274),seq(277,291),seq(296,299),seq(301,304),seq(306,309),seq(313,315),seq(318,320),321,355,363)
#howland:iseq <- c(seq(182,193),seq(195,203),seq(206,208),211,213,seq(215,217),224,230,231,seq(233,236),seq(244,254),seq(256,260),seq(262,258),seq(271,274),seq(277,291),seq(296,299),seq(301,304),seq(313,315),318,355,363)
#iseq <- c(seq(183,187),189,190,191,192,196,197,198,199,200,202,207,210,211,212,213,214,218,222,225,228,229,233,234,235,236,237,238,239,240,242,243,244,245,247,250,251,252,253,254,255,256,259,260,seq(266,272),seq(274,278),280,283,285,286,287,290,291,292,293,294,295,300,301,303,304,308,311,312,313,314,315,316,319,321,355)

#iseq <- as.character(c(seq(1,333,1),seq(348,364,1)))
#i <- "040"
#i <- "028"
#i <- "131"
iseq <- as.character(c(1,10,12,15,17,20,21,22,46,57,93,125,122,131,137,150,148))
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
  if(as.numeric(i)<182){
    year <- 2018
  }
  else{
    year <- 2017
  }
  fileName <- paste("dailyNDVI_GOES/","GOES_Diurnal_",siteName,"_",year,i,".csv",sep="")
  print(fileName)
  dat <- read.csv(fileName,header=FALSE)
  data <- list()
  print(dim(dat))
  data$x <- as.numeric(dat[3,])
  data$y <- as.numeric(dat[2,])
  if(TZ==6){
    solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(i),origin=as.Date(paste((year-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/Chicago"))$solarNoon
    solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
  }
  data$k <- solarNoonTime
  plot(data$x,data$y,pch=20,ylim=c(0,1))
  xseq <- seq(5,20,0.01)
  #points(xseq,diurnalExp(a=0.004,c=0.2,k=12,xseq=xseq),col="red")
  
  outFileName <- paste(siteName,"_",as.character(i),"_varBurn5.RData",sep="")
  if(!file.exists(outFileName)){
    j.model <- createBayesModel.Diurnal(siteName=siteName,data)
    var.burn <- runMCMC_Model(j.model = j.model,variableNames=c("a","c","prec","p.cloud","alp","bet"),baseNum = 500000,iterSize=100000,maxGBR=20,maxIter=10000000)#,baseNum = 1000000,iterSize = 70000)
    if(typeof(var.burn)!=typeof(FALSE)){
      save(var.burn,file=outFileName)
    }
  }
}
