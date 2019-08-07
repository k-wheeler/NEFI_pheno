###Spring Forecast Uncertainty/Sensitivity Analysis
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary") 
library(rjags)
library(runjags)
library("PhenoForecast")
library("PhenologyBayesModeling")
library(ecoforecastR)

siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
forecastDataFolder <- "PhenologyForecastData/ForecastOutputs/AllForecasts/"
dataDirectory <- "PhenologyForecastData/"
Nmc <- 1000 #Number of model runs
#allDates <- seq(as.Date("2019-01-23"),as.Date("2019-06-06"),"day")
allDates <- c(seq(as.Date("2019-01-23"),as.Date("2019-01-25"),"day"),
              as.Date("2019-02-03"),as.Date("2019-02-05"),
              seq(as.Date("2019-02-07"),as.Date("2019-06-06"),"day"))
allDates <- allDates[1:(length(allDates)/2)]
#allDates <- allDates[(length(allDates)/2):length(allDates)]
i <- 1
#allDates <- allDates[(length(allDates)/2):length(allDates)]
##General site-specific info
siteName <- as.character(siteData$siteName[i])
print(siteName)
URL <- as.character(siteData$URL[i])
lat <- as.numeric(siteData$Lat[i])
long <- as.numeric(siteData$Long[i])
station <- as.character(siteData$metStation[i])
siteStartDate <- as.character(siteData$startDate[i])
pdf(file=paste(siteName,"_PhenologyForecast_uncertaintyAnalysis.pdf",sep=""),height=10,width=6)
par(mfrow=c(3,1))

##Date info 
for(d in 1:length(allDates)){
  
  calEndDate <- allDates[d]
  #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
  #calEndDate <- as.Date("2019-05-20")
  print(calEndDate)
  forStartDate <- calEndDate + 1
  forEndDate <- (forStartDate+13)
  forDates <- seq(forStartDate,forEndDate,"day")
  yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
  
  forecastLength <- 14
  plotDates <- as.numeric(format(forDates,"%j"))
  
  ##Load data for plots
  URL <- as.character(siteData$URL[i])
  URL2 <- as.character(siteData$URL2[i])
  URL3 <- as.character(siteData$URL3[i])
  if(nchar(URL2)>0){
    URL <- c(URL,URL2)
    if(nchar(URL3)>0){
      URL <- c(URL,URL3)
    }
  }
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
  startDate <- as.Date(siteData[i,7])
  rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_PC.csv",sep="")
  rescaleData <- read.csv(rescaleFile,header=TRUE)
  cMeans.p <- rescaleData$cMeans.p
  dMeans.p <- rescaleData$dMeans.p
  rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_MN.csv",sep="")
  rescaleData <- read.csv(rescaleFile,header=TRUE)
  cMeans.mn <- rescaleData$cMeans.mn
  dMeans.mn <- rescaleData$dMeans.mn
  rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_ME.csv",sep="")
  rescaleData <- read.csv(rescaleFile,header=TRUE)
  cMeans.me <- rescaleData$cMeans.me
  dMeans.me <- rescaleData$dMeans.me
  
  phenoData <- matrix(nrow=0,ncol=32)
  for(u in 1:length(URL)){
    phenoDataSub <- download.phenocam(URL[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }
  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  phenoData <- phenoData[phenoData$date<calEndDate,]
  
  p <- phenoData$gcc_mean[phenoData$year==2019]
  
  time.p.date <- as.Date(phenoData$date[phenoData$year==2019])
  time.p <-  as.numeric(format(time.p.date,"%j"))
  
  p <- rescale(c=cMeans.p[length(cMeans.p)],d=dMeans.p[length(dMeans.p)],yseq=p)
  p[p<0] <- 0
  
  ##MODIS NDVI and EVI
  mn <- prepareMODIS(startDate=startDate,endDate=calEndDate,metric="NDVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
  me <- prepareMODIS(startDate=startDate,endDate=calEndDate,metric="EVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
  mn <- rescale(yseq=mn,c=cMeans.mn[length(cMeans.mn)],d=dMeans.mn[length(dMeans.mn)])
  me <- rescale(yseq=me,c=cMeans.me[length(cMeans.me)],d=dMeans.me[length(dMeans.me)])
  mn[mn<0] <- 0
  me[me<0] <- 0
  
  ############Random Walk
  RWfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_randomWalk_outBurn.RData",sep="")
  #RWfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_randomWalk_outBurn.RData",sep="")
  if(file.exists(RWfileName)){
    load(RWfileName)
    out.mat.par <- data.frame(as.matrix(outBurnRW$params))
    dayNumber <- dim(as.matrix(outBurnRW$predict))[2]-14
    out.mat.RW <- as.matrix(outBurnRW$predict)
    out.mat.RW <- out.mat.RW[,1:dayNumber]
    ICs <- out.mat.RW[,dayNumber]
    
    ##Deterministic prediction
    RW.det <- forecastRW(IC=mean(ICs),
                         Q=0,
                         n=1,
                         NT=14)
    
    ##Initial Condition Uncertainty
    rndNums <- sample.int(nrow(out.mat.par),Nmc,replace=TRUE)
    RW.IC <- forecastRW(IC=ICs[rndNums],
                        Q=0,
                        n=Nmc,
                        NT=14)
    RW.IC.ci <- apply(RW.IC,2,quantile,c(0.025,0.5,0.975))
    
    ##Process Uncertainty
    Qmc <- 1/sqrt(out.mat.par[rndNums,"p.proc"])
    RW.IPDE <- forecastRW(IC=ICs[rndNums],
                          Q=Qmc,
                          n=Nmc,
                          NT=14)
    RW.IPDE.ci <- apply(RW.IPDE,2,quantile,c(0.025,0.5,0.975))
    plotRun(out.mat = out.mat.RW,forecastType = "randomWalk",endDate=calEndDate)
    ecoforecastR::ciEnvelope(plotDates,RW.IPDE.ci[1,],RW.IPDE.ci[3,],col=adjustcolor("blue",0.8))
    ecoforecastR::ciEnvelope(plotDates,RW.IC.ci[1,],RW.IC.ci[3,],col=adjustcolor("gray",1))
    lines(plotDates,RW.det,col="purple",lwd=2)
    abline(v=as.numeric(format(calEndDate,"%j")),col="red",lwd=2)
  
    points(time.p,p,pch=20,col="red")
    points(time.p,mn,col="blue",pch=3,cex=2)
    points(time.p,me,col="blue",pch=1,cex=2)
    }
  #################Basic Logistic
  logFileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_logistic_outBurn.RData",sep="")
  #logFileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_logistic_outBurn.RData",sep="")
  if(file.exists(logFileName)){
    load(logFileName)
    out.mat.par <- data.frame(as.matrix(outBurnL$params))
    dayNumber <- dim(as.matrix(outBurnL$predict))[2]-14
    out.mat.L <- as.matrix(outBurnL$predict)
    out.mat.L <- out.mat.L[,1:dayNumber]
    ICs <- out.mat.L[,dayNumber]
    
    ##Deterministic prediction
    L.det <- forecastLog(IC=mean(ICs),
                         r=mean(out.mat.par$r),
                         Q=0,
                         n=1,
                         NT=14)
    
    ##Initial Condition Uncertainty
    rndNums <- sample.int(nrow(out.mat.par),Nmc,replace=TRUE)
    L.I <- forecastLog(IC=ICs[rndNums],
                       r=mean(out.mat.par$r),
                       Q=0,
                       n=Nmc,
                       NT=14)
    L.I.ci <- apply(L.I,2,quantile,c(0.025,0.5,0.975))
    
    ##Parameter Uncertainty <- Don't have for RW (will have for log and logCov)
    L.IP <- forecastLog(IC=ICs[rndNums],
                        r=out.mat.par$r[rndNums],
                        Q=0,
                        n=Nmc,
                        NT=14)
    L.IP.ci = apply(L.IP,2,quantile,c(0.025,0.5,0.975))
    
    ##Process Uncertainty
    Qmc <- 1/sqrt(out.mat.par[rndNums,"p.proc"])
    L.IPDE <- forecastLog(IC=ICs[rndNums],
                          r=out.mat.par$r[rndNums],
                          Q=Qmc,
                          n=Nmc,
                          NT=14)
    L.IPDE.ci <- apply(L.IPDE,2,quantile,c(0.025,0.5,0.975))
    
    plotRun(out.mat = out.mat.L,forecastType = "logistic",endDate=calEndDate)
    ecoforecastR::ciEnvelope(plotDates,L.IPDE.ci[1,],L.IPDE.ci[3,],col=adjustcolor("blue",0.8))
    ecoforecastR::ciEnvelope(plotDates,L.IP.ci[1,],L.IP.ci[3,],col=adjustcolor("green",0.8))
    ecoforecastR::ciEnvelope(plotDates,L.I.ci[1,],L.I.ci[3,],col=adjustcolor("gray",1))
    lines(plotDates,L.det,col="purple",lwd=2)
    abline(v=as.numeric(format(calEndDate,"%j")),col="red",lwd=3)
  
    points(time.p,p,pch=20,col="red")
    points(time.p,mn,col="blue",pch=3,cex=2)
    points(time.p,me,col="blue",pch=1,cex=2)
    }
  
  ###################Logistic with Covariates
  LCfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
  if(file.exists(LCfileName)){
    ##Assemble Sf forecast
    GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/",siteName,"/",calEndDate,"/",sep="")
    GEFS_Files <- dir(path=GEFS_Directory,pattern="NOAA_GEFS")
    
    SfsALL <- createSf(lat=lat,long=long,dates=yrDates,siteName=siteName,
             dataDirectory=dataDirectory,endDate=forEndDate,
             GEFS_Files=GEFS_Files,GEFS_Directory=GEFS_Directory,
             forecastLength=forecastLength, station=station,
             calDatesT = FALSE) 
    
    #save(SfsALL,file="SfsALL.RData")
    #load("SfsALL_140.RData")
    #load("SfsALL_140_new.RData")
    #load("harvard_2009-01-01_2019-05-20_LC3_outBurn.RData")
    #outBurnLC2 <- outBurnLC
    #Sf.means <- matrix(colMeans(SfsALL),1,forecastLength)
    
    #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
    
    load(LCfileName)
    out.mat.par <- data.frame(as.matrix(outBurnLC2$params))
    dayNumber <- dim(as.matrix(outBurnLC2$predict))[2]-14
    out.mat.LC <- as.matrix(outBurnLC2$predict)
    out.mat.LC <- out.mat.LC[,1:dayNumber]
    ICs <- out.mat.LC[,dayNumber]
    
    ##Deterministic prediction
    LC.det <- forecastLogCov(IC=mean(ICs),
                             trans = mean(out.mat.par$trans),
                             b1 = mean(out.mat.par$b1),
                             Sf = matrix(SfsALL$Sf[(length(SfsALL$Sf)-13):length(SfsALL$Sf)],1,14),
                             Q=0,
                             n=1,
                             NT=14)
    
    ##Initial Condition Uncertainty
    rndNums <- sample.int(nrow(out.mat.par),Nmc,replace=TRUE)
    LC.I <- forecastLogCov(IC=ICs[rndNums],
                           trans = mean(out.mat.par$trans),
                           b1 = mean(out.mat.par$b1),
                           Sf = matrix(SfsALL$Sf[(length(SfsALL$Sf)-13):length(SfsALL$Sf)],1,14),
                           Q=0,
                           n=Nmc,
                           NT=14)
    LC.I.ci <- apply(LC.I,2,quantile,c(0.025,0.5,0.975))
    
    ##Parameter Uncertainty 
    LC.IP <- forecastLogCov(IC=ICs[rndNums],
                            trans = out.mat.par$trans[rndNums],
                            b1 = out.mat.par$b1[rndNums],
                            Sf = matrix(SfsALL$Sf[(length(SfsALL$Sf)-13):length(SfsALL$Sf)],1,14),
                            Q=0,
                            n=Nmc,
                            NT=14)
    LC.IP.ci = apply(LC.IP,2,quantile,c(0.025,0.5,0.975))
    
    ##Driver Uncertainty
    SfsSamp <- matrix(nrow=Nmc,ncol=0)
    #SfsALL
    for(s in (length(SfsALL$Sf)-13):length(SfsALL$Sf)){
      SfsSamp <- cbind(SfsSamp,rnorm(Nmc,SfsALL$Sf[s],1/sqrt(SfsALL$Sfprec[s])))
    }
    LC.IPD <- forecastLogCov(IC=ICs[rndNums],
                             trans = out.mat.par$trans[rndNums],
                             b1 = out.mat.par$b1[rndNums],
                             Sf = SfsSamp,
                             Q=0,
                             n=Nmc,
                             NT=14)
    LC.IPD.ci = apply(LC.IPD,2,quantile,c(0.025,0.5,0.975))
    
    ##Process Uncertainty
    Qmc <- 1/sqrt(out.mat.par[rndNums,"p.proc"])
    LC.IPDE <- forecastLogCov(IC=ICs[rndNums],
                              trans = out.mat.par$trans[rndNums],
                              b1 = out.mat.par$b1[rndNums],
                              Sf = SfsSamp,
                              Q=Qmc,
                              n=Nmc,
                              NT=14)
    LC.IPDE.ci <- apply(LC.IPDE,2,quantile,c(0.025,0.5,0.975))
    
    plotRun(out.mat = out.mat.LC,forecastType = "logistic covariate",endDate=calEndDate)
    ecoforecastR::ciEnvelope(plotDates,LC.IPDE.ci[1,],LC.IPDE.ci[3,],col=adjustcolor("blue",0.8))
    ecoforecastR::ciEnvelope(plotDates,LC.IPD.ci[1,],LC.IPD.ci[3,],col=adjustcolor("pink",0.8))
    ecoforecastR::ciEnvelope(plotDates,LC.IP.ci[1,],LC.IP.ci[3,],col=adjustcolor("green",0.8))
    ecoforecastR::ciEnvelope(plotDates,LC.I.ci[1,],LC.I.ci[3,],col=adjustcolor("gray",1))
    lines(plotDates,LC.det,col="purple",lwd=2)
    abline(v=as.numeric(format(calEndDate,"%j")),col="red",lwd=3)
  
    points(time.p,p,pch=20,col="red")
    points(time.p,mn,col="blue",pch=3,cex=2)
    points(time.p,me,col="blue",pch=1,cex=2)
    }
}

dev.off()


