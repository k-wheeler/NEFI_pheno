###Spring Forecast Uncertainty/Sensitivity Analysis
install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary") 
library(rjags)
library(runjags)
library("PhenoForecast")
library("PhenologyBayesModeling")
library(ecoforecastR)

siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
forecastDataFolder <- "PhenologyForecastData/ForecastOutputs/AllForecasts/"
Nmc <- 1000 #Number of model runs
#allDates <- seq(as.Date("2019-01-23"),as.Date("2019-06-06"),"day")
allDates <- c(seq(as.Date("2019-01-23"),as.Date("2019-01-25"),"day"),
              as.Date("2019-02-03"),as.Date("2019-02-05"),
              seq(as.Date("2019-02-07"),as.Date("2019-06-06"),"day"))
i <- 1

##General site-specific info
siteName <- as.character(siteData$siteName[i])
print(siteName)
URL <- as.character(siteData$URL[i])
lat <- as.numeric(siteData$Lat[i])
long <- as.numeric(siteData$Long[i])
station <- as.character(siteData$metStation[i])
siteStartDate <- as.character(siteData$startDate[i])
pdf(file=paste(siteName,"_PhenologyForecast_uncertaintyAnalysis.pdf",sep=""),height=6,width=10)
par(mfrow=c(4,2))

##Date info 
for(d in 1:length(allDates)){
  
  calEndDate <- allDates[d]
  print(calEndDate)
  forStartDate <- calEndDate + 1
  forEndDate <- (forStartDate+13)
  forDates <- seq(forStartDate,forEndDate,"day")
  
  forecastLength <- 14
  plotDates <- as.numeric(format(forDates,"%j"))
  
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
    ecoforecastR::ciEnvelope(plotDates,RW.IC.ci[1,],RW.IC.ci[3,],col=adjustcolor("gray",0.8))
    lines(plotDates,RW.det,col="purple",lwd=2)
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
    ecoforecastR::ciEnvelope(plotDates,L.I.ci[1,],L.I.ci[3,],col=adjustcolor("gray",0.8))
    lines(plotDates,L.det,col="purple",lwd=2)
  }
  
  ###################Logistic with Covariates
  LCfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
  if(file.exists(LCfileName)){
    ##Assemble Sf forecast
    GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/",siteName,"/",calEndDate,"/",sep="")
    GEFS_Files <- dir(path=GEFS_Directory,pattern="NOAA_GEFS")
    
    TairsForecast <- numeric()
    for(e in 1:length(GEFS_Files)){
      TairsForecastInd <- load_GEFS_Forecast(dataDirectory=GEFS_Directory,fileName=GEFS_Files[e])
      TairsForecast <- cbind(TairsForecast,TairsForecastInd)
    }
    SfsALL <- matrix(nrow=0,ncol=length(forDates))
    
    for(e in 1:ncol(TairsForecast)){
      Sfs <- calSf(Tairs=TairsForecast[,e],dates=forDates)
      SfsALL <- rbind(SfsALL,Sfs)
    }
    
    #save(SfsALL,file="SfsALL.RData")
    #load("SfsALL.RData")
    
    Sf.means <- matrix(colMeans(SfsALL),1,forecastLength)
    
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
                             Sf = Sf.means,
                             Q=0,
                             n=1,
                             NT=14)
    
    ##Initial Condition Uncertainty
    rndNums <- sample.int(nrow(out.mat.par),Nmc,replace=TRUE)
    LC.I <- forecastLogCov(IC=ICs[rndNums],
                           trans = mean(out.mat.par$trans),
                           b1 = mean(out.mat.par$b1),
                           Sf = Sf.means,
                           Q=0,
                           n=Nmc,
                           NT=14)
    LC.I.ci <- apply(LC.I,2,quantile,c(0.025,0.5,0.975))
    
    ##Parameter Uncertainty <- Don't have for RW (will have for log and logCov)
    LC.IP <- forecastLogCov(IC=ICs[rndNums],
                            trans = out.mat.par$trans[rndNums],
                            b1 = out.mat.par$b1[rndNums],
                            Sf = Sf.means,
                            Q=0,
                            n=Nmc,
                            NT=14)
    LC.IP.ci = apply(LC.IP,2,quantile,c(0.025,0.5,0.975))
    
    ##Driver Uncertainty <- Don't have for RW (will have for logCov)
    drow = sample.int(nrow(SfsALL),Nmc,replace=TRUE)
    LC.IPD <- forecastLogCov(IC=ICs[rndNums],
                             trans = out.mat.par$trans[rndNums],
                             b1 = out.mat.par$b1[rndNums],
                             Sf = SfsALL[drow,],
                             Q=0,
                             n=Nmc,
                             NT=14)
    LC.IPD.ci = apply(LC.IPD,2,quantile,c(0.025,0.5,0.975))
    
    ##Process Uncertainty
    Qmc <- 1/sqrt(out.mat.par[rndNums,"p.proc"])
    LC.IPDE <- forecastLogCov(IC=ICs[rndNums],
                              trans = out.mat.par$trans[rndNums],
                              b1 = out.mat.par$b1[rndNums],
                              Sf = SfsALL[drow,],
                              Q=Qmc,
                              n=Nmc,
                              NT=14)
    LC.IPDE.ci <- apply(LC.IPDE,2,quantile,c(0.025,0.5,0.975))
    plotRun(out.mat = out.mat.LC,forecastType = "logistic",endDate=calEndDate)
    ecoforecastR::ciEnvelope(plotDates,LC.IPDE.ci[1,],LC.IPDE.ci[3,],col=adjustcolor("blue",0.8))
    ecoforecastR::ciEnvelope(plotDates,LC.IPDE.ci[1,],LC.IPD.ci[3,],col=adjustcolor("pink",0.8))
    ecoforecastR::ciEnvelope(plotDates,LC.IP.ci[1,],LC.IP.ci[3,],col=adjustcolor("green",0.8))
    ecoforecastR::ciEnvelope(plotDates,LC.I.ci[1,],LC.I.ci[3,],col=adjustcolor("gray",0.8))
    lines(plotDates,LC.det,col="purple",lwd=2)
  }
}

dev.off()
