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
Nmc <- 10000 #Number of model runs
#allDates <- seq(as.Date("2019-01-23"),as.Date("2019-06-06"),"day")
allDates <- c(seq(as.Date("2019-01-23"),as.Date("2019-01-25"),"day"),
              as.Date("2019-02-03"),as.Date("2019-02-05"),
              seq(as.Date("2019-02-07"),as.Date("2019-06-06"),"day"))
#allDates <- allDates[1:(length(allDates)/2)]
#allDates <- allDates[(length(allDates)/2):length(allDates)]
i <- 6

##General site-specific info
siteName <- as.character(siteData$siteName[i])
print(siteName)
plotFolder <- paste("VPplots/",siteName,"/",sep="")
URL <- as.character(siteData$URL[i])
lat <- as.numeric(siteData$Lat[i])
long <- as.numeric(siteData$Long[i])
station <- as.character(siteData$metStation[i])
siteStartDate <- as.character(siteData$startDate[i])
#pdf(file=paste(siteName,"_PhenologyForecast_variancePartition.pdf",sep=""),height=10,width=6)
#par(mfrow=c(3,1))

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

    varI     <- apply(RW.IC,2,var)
    varIPDE  <- apply(RW.IPDE,2,var)

    varMat   <- rbind(varI,varIPDE)
    
    ## out-of-sample stacked area plot
    V.pred.rel <- apply(varMat,2,function(x) {x/max(x)})
    plotFileName <- paste(plotFolder,"randomWalk/",siteName,"_randomWalk_",calEndDate,"_varPartition.png",sep="")
    png(file=plotFileName, width=10, height=5,units="in",res=1000)
    par(mai=c(1,1.2,0.5,0.5))
    plot(forDates,V.pred.rel[1,],ylim=c(0,1),type='n',main="",ylab="Proportion of Variance",xlab="Time",cex.lab=2,cex.axis=2)
    ciEnvelope(forDates,rep(0,ncol(V.pred.rel)),V.pred.rel[1,],col="blue")
    ciEnvelope(forDates,V.pred.rel[1,],V.pred.rel[2,],col="gray")
    legend("topleft",legend=calEndDate,pch=20,col="white",cex=1.5)
    dev.off()

    #legend("topleft",legend=c("RandomEffect","Process","Driver","Parameter","InitCond"),col=rev(N.cols),lty=1,lwd=5)
    
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
    
    varI     <- apply(L.I,2,var)
    varIP    <- apply(L.IP,2,var)
    varIPDE  <- apply(L.IPDE,2,var)
    varMat   <- rbind(varI,varIP,varIPDE)
    
    ## out-of-sample stacked area plot
    V.pred.rel <- apply(varMat,2,function(x) {x/max(x)})
    plotFileName <- paste(plotFolder,"logistic/",siteName,"_logistic_",calEndDate,"_varPartition.png",sep="")
    png(file=plotFileName, width=10, height=5,units="in",res=1000)
    par(mai=c(1,1.2,0.5,0.5))
    plot(forDates,V.pred.rel[1,],ylim=c(0,1),type='n',main="",ylab="Proportion of Variance",xlab="Time",cex.lab=2,cex.axis=2)
    ciEnvelope(forDates,rep(0,ncol(V.pred.rel)),V.pred.rel[1,],col="gray")
    ciEnvelope(forDates,V.pred.rel[1,],V.pred.rel[2,],col="green")
    ciEnvelope(forDates,V.pred.rel[2,],V.pred.rel[3,],col="blue")
    legend("topleft",legend=calEndDate,pch=20,col="white",cex=1.5)
    dev.off()

    #legend("topleft",legend=c("RandomEffect","Process","Driver","Parameter","InitCond"),col=rev(N.cols),lty=1,lwd=5)
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
    
    save(SfsALL,file=paste("SfsALL_",calEndDate,".RData",sep=""))
    #load("SfsALL_140.RData")
    #load("SfsALL_140_new.RData")
    #load("harvard_2009-01-01_2019-05-20_LC2_outBurn.RData")
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
    for(s in (length(SfsALL$Sf)-13):length(SfsALL$Sf)){ ##Only need the last 14
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
    
    varI     <- apply(LC.I,2,var)
    varIP    <- apply(LC.IP,2,var)
    varIPD   <- apply(LC.IPD,2,var)
    varIPDE  <- apply(LC.IPDE,2,var)
    varMat   <- rbind(varI,varIP,varIPD,varIPDE)
    
    ## out-of-sample stacked area plot
    V.pred.rel <- apply(varMat,2,function(x) {x/max(x)})
    plotFileName <- paste(plotFolder,"LC2/",siteName,"_LC2_",calEndDate,"_varPartition.png",sep="")
    png(file=plotFileName, width=10, height=5,units="in",res=1000)
    par(mai=c(1,1.2,0.5,0.5))
    plot(forDates,V.pred.rel[1,],ylim=c(0,1),type='n',main="",ylab="Proportion of Variance",xlab="Time",cex.lab=2,cex.axis=2)
    ciEnvelope(forDates,V.pred.rel[3,],V.pred.rel[4,],col="blue")
    ciEnvelope(forDates,rep(0,ncol(V.pred.rel)),V.pred.rel[1,],col="gray")
    ciEnvelope(forDates,V.pred.rel[1,],V.pred.rel[2,],col="green")
    ciEnvelope(forDates,V.pred.rel[2,],V.pred.rel[3,],col="pink")
    legend("topleft",legend=calEndDate,pch=20,col="white",cex=1.5)
    dev.off()
    #legend("topleft",legend=c("RandomEffect","Process","Driver","Parameter","InitCond"),col=rev(N.cols),lty=1,lwd=5)
  }
}

dev.off()
