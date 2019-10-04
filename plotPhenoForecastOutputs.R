install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary") 

library("PhenoForecast")
library("PhenologyBayesModeling")
library("coda")
library("dplyr")
library("rjags")

##Read in data
#siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
#dataDirectory="PhenologyForecastData/"
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 15

endDate <- (Sys.Date()-2)
#startDate <- as.Date("2013-01-01")
#endDate <- as.Date("2019-01-27")
i <- 10
i <- 1
iseq <- c(seq(1,4),6,8,9,10,15,16,seq(18,27))
iseq <- c(1)
iseq <- c(1,2,3,4,6,15,16,20,18,24)
outFileName <- paste("PhenologyForecastOutput_allSites_",endDate,".pdf",sep="")
#iseq <- c(seq(1,6),8,9,10,15,16,seq(18,27))
pdf(outFileName,height=6,width=10)

for(i in iseq){  
  siteName <- as.character(siteData[i,1])
  print(siteName)
  saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
  URL <- as.character(siteData$URL[i])
  URL2 <- as.character(siteData$URL2[i])
  URL3 <- as.character(siteData$URL3[i])
  if(nchar(URL2)>0){
    URL <- c(URL,URL2)
    if(nchar(URL3)>0){
      URL <- c(URL,URL3)
    }
  }
  #print(URL)
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
  startDate <- as.Date(siteData[i,7])
  
  ##Load rescale data
  rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_PC.csv",sep="")
  if(file.exists(rescaleFile)){
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
    
    ##Random Walk:
    RWFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_RW_outBurn.RData",sep="")
    if(file.exists(RWFile)){
      load(RWFile)
      
      ##Plot Parameter Density Plots
      out.mat.par <- data.frame(as.matrix(outBurnRW$params))
      print(colnames(out.mat.par))
      par(mfrow=c(2,2))
      plot(density(out.mat.par$p.ME),main="Density of p.ME")
      plot(density(out.mat.par$p.MN),main="Density of p.MN")
      plot(density(out.mat.par$p.PC),main="Density of p.PC")
      plot(density(out.mat.par$p.proc),main="Density of p.proc")
      
      plot(density(sqrt(1/out.mat.par$p.PC)),main="Density of PC Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.ME)),main="Density of ME Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.MN)),main="Density of MN Ob. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.proc)),main="Density of Process Error (SD)")
      
      par(mfrow=c(1,1))
      ##Plot Forecast including calibration
      dayNumber <- dim(as.matrix(outBurnRW$predict))[2]
      out.mat.RW <- as.matrix(outBurnRW$predict)
      print(length(seq(1,dayNumber,1)))
      print(dim(out.mat.RW))
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                         out.mat=out.mat.RW,forecastType = "randomWalk",days=seq(1,dayNumber,1),
                         endDate = endDate)#,dates=seq(startDate,(endDate+forecastLength),"day"))
      offset <- 365-as.numeric(format(startDate,"%j"))
      for(i in seq(offset,dayNumber,365)){
        abline(v=i,col="red")
      }
      abline(v=(dayNumber-forecastLength+1),col="purple")
      print("Done all")
      
      ##Plot Current year
      lengthLastYear <- (as.numeric(format(endDate,"%j"))+forecastLength)
      print(lengthLastYear)
      lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
      out.mat.lastYear <- out.mat.RW[,lastYearIndices]
      
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                         out.mat=out.mat.lastYear,forecastType = "randomWalk",
                         days=seq(1,lengthLastYear,1),endDate=endDate)#,
      #dates=seq(as.Date("2019-07-01"),(endDate+forecastLength),"day"))
      abline(v=(lengthLastYear-forecastLength+1),col="purple")
      
      ##Add on data:
      ##PhenoCam
      phenoData <- matrix(nrow=0,ncol=32)
      for(u in 1:length(URL)){
        print(URL[u])
        phenoDataSub <- download.phenocam(URL[u])
        phenoData <- rbind(phenoData,phenoDataSub)
      }
      ##Order and remove duplicate PC data
      phenoData2 <- phenoData[order(phenoData$date),]
      phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
      phenoData <- phenoData3
      phenoData <- phenoData[phenoData$date<endDate,]
      
      p <- phenoData$gcc_mean[phenoData$year==2019]
      
      time.p.date <- as.Date(phenoData$date[phenoData$year==2019])
      time.p <-  as.numeric(format(time.p.date,"%j"))
      
      p <- rescale(c=cMeans.p[length(cMeans.p)],d=dMeans.p[length(dMeans.p)],yseq=p)
      points(time.p,p,pch=20,col="red")
      
      ##MODIS NDVI and EVI
      mn <- prepareMODIS(startDate=startDate,endDate=endDate,metric="NDVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
      me <- prepareMODIS(startDate=startDate,endDate=endDate,metric="EVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
      mn <- rescale(yseq=mn,c=cMeans.mn[length(cMeans.mn)],d=dMeans.mn[length(dMeans.mn)])
      me <- rescale(yseq=me,c=cMeans.me[length(cMeans.me)],d=dMeans.me[length(dMeans.me)])
      points(time.p,mn,col="green",pch="+")
      points(time.p,me,col="green",pch=1)
      legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
    }
    ##Basic logistic:
    LFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_log_outBurn.RData",sep="")
    if(file.exists(LFile)){
      load(LFile)
      out.mat.par <- data.frame(as.matrix(outBurnL$params))
      print(colnames(out.mat.par))
      par(mfrow=c(2,3))
      plot(density(out.mat.par$p.ME),main="Density of p.ME")
      plot(density(out.mat.par$p.MN),main="Density of p.MN")
      plot(density(out.mat.par$p.PC),main="Density of p.PC")
      plot(density(out.mat.par$p.proc),main="Density of p.proc")
      plot(density(out.mat.par$r),main="Density of r")
      par(mfrow=c(2,3))
      plot(density(sqrt(1/out.mat.par$p.PC)),main="Density of PC Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.ME)),main="Density of ME Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.MN)),main="Density of MN Ob. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.proc)),main="Density of Process Error (SD)")
      plot(density(out.mat.par$r),main="Density of r")
      
      ##Plot
      par(mfrow=c(1,1))
      dayNumber <- dim(as.matrix(outBurnL$predict))[2]
      out.mat.L <- as.matrix(outBurnL$predict)
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                         out.mat=out.mat.L,forecastType = "Logistic",days=seq(1,dayNumber,1),
                         endDate = endDate)#,dates=seq(startDate,(endDate+forecastLength),"day"))#2484
      
      for(i in seq(1,2484,182)){
        abline(v=i,col="red")
      }
      abline(v=(dayNumber-forecastLength+1),col="purple")
      #points(time.p,p,pch=20,col="red")
      #points(time.p,mn,col="green",pch=3)
      #points(time.p,me,col="green",pch=1)
      legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
      
      
      ##Plot Current year
      lengthLastYear <- lengthLastYear - 212 ##Done for autumn
      lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
      out.mat.lastYear <- out.mat.L[,lastYearIndices]
      
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                         out.mat=out.mat.lastYear,forecastType = "Logistic",days=seq(213,(lengthLastYear+212),1),
                         endDate=endDate,xlim=c(215,(lengthLastYear+220)))#,
      #dates=seq(as.Date("2019-07-01"),(endDate+forecastLength),"day"))
      
      abline(v=(lengthLastYear-forecastLength+1),col="purple")
      
      ##Add on data points
      points(time.p,p,pch=20,col="red")
      points(time.p,mn,col="green",pch=3)
      points(time.p,me,col="green",pch=1)
      legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
    }
    ##Logistic with covariates model
    LCFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_20_meanTair_outBurn.RData",sep="")
    print(LCFile)
    if(file.exists(LCFile)){
      print("file exists")
      load(LCFile)
      outBurnLC <- outBurn
      out.mat.par <- data.frame(as.matrix(outBurnLC$params))
      print(colnames(out.mat.par))
      par(mfrow=c(2,3))
      plot(density(out.mat.par$p.ME),main="Density of p.ME")
      plot(density(out.mat.par$p.MN),main="Density of p.MN")
      plot(density(out.mat.par$p.PC),main="Density of p.PC")
      plot(density(out.mat.par$p.proc),main="Density of p.proc")
      plot(density(out.mat.par$MOF),main="Density of MOF")
      plot(density(out.mat.par$fallLength),main="Density of fallLength")
      plot(density(out.mat.par$sSlope),main="Density of sSlope")
      
      par(mfrow=c(2,3))
      plot(density(sqrt(1/out.mat.par$p.PC)),main="Density of PC Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.ME)),main="Density of ME Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.MN)),main="Density of MN Ob. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.proc)),main="Density of Process Error (SD)")
      
      plot(density(out.mat.par$b1),main="Density of b1")
      #plot(density(out.mat.par$b0),main="Density of b0")
      
      ##Plot
      par(mfrow=c(1,1))
      dayNumber <- dim(as.matrix(outBurnLC$predict))[2]
      out.mat.LC <- as.matrix(outBurnLC$predict)
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                         out.mat=out.mat.LC,forecastType = "CDD_20",days=seq(1,dayNumber,1),
                         endDate = endDate)#2484
      
      for(i in seq(1,2484,182)){
        abline(v=i,col="red")
      }
      abline(v=(dayNumber-forecastLength+1),col="purple")
      
      ##Plot Current year
      lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
      out.mat.lastYear <- out.mat.LC[,lastYearIndices]
      
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.lastYear,
                         forecastType = "CDD_20",days=seq(213,(lengthLastYear+212),1),
                         endDate = endDate,xlim=c(215,(lengthLastYear+220)))
      abline(v=(lengthLastYear-forecastLength+1),col="purple")
      
      ##Add on data points
      points(time.p,p,pch=20,col="red")
      points(time.p,mn,col="green",pch=3)
      points(time.p,me,col="green",pch=1)
      legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
    }
    ##Logistic with covariates model
    LCFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_NA_meanTair_outBurn.RData",sep="")
    if(file.exists(LCFile)){
      load(LCFile)
      outBurnLC <- outBurn
      out.mat.par <- data.frame(as.matrix(outBurnLC$params))
      print(colnames(out.mat.par))
      par(mfrow=c(2,3))
      plot(density(out.mat.par$p.ME),main="Density of p.ME")
      plot(density(out.mat.par$p.MN),main="Density of p.MN")
      plot(density(out.mat.par$p.PC),main="Density of p.PC")
      plot(density(out.mat.par$p.proc),main="Density of p.proc")
      plot(density(out.mat.par$MOF),main="Density of MOF")
      plot(density(out.mat.par$fallLength),main="Density of fallLength")
      plot(density(out.mat.par$sSlope),main="Density of sSlope")
      plot(density(out.mat.par$baseTemp),main="Density of baseTemp")
      
      par(mfrow=c(2,3))
      plot(density(sqrt(1/out.mat.par$p.PC)),main="Density of PC Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.ME)),main="Density of ME Obs. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.MN)),main="Density of MN Ob. Error (SD)")
      plot(density(sqrt(1/out.mat.par$p.proc)),main="Density of Process Error (SD)")
      
      plot(density(out.mat.par$b1),main="Density of b1")
      #plot(density(out.mat.par$b0),main="Density of b0")
      
      ##Plot
      par(mfrow=c(1,1))
      dayNumber <- dim(as.matrix(outBurnLC$predict))[2]
      out.mat.LC <- as.matrix(outBurnLC$predict)
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                         out.mat=out.mat.LC,forecastType = "CDD_NA",days=seq(1,dayNumber,1),
                         endDate = endDate)#2484
      
      for(i in seq(1,2484,182)){
        abline(v=i,col="red")
      }
      abline(v=(dayNumber-forecastLength+1),col="purple")
      
      ##Plot Current year
      lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
      out.mat.lastYear <- out.mat.LC[,lastYearIndices]
      
      plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.lastYear,
                         forecastType = "CDD_NA",days=seq(213,(lengthLastYear+212),1),
                         endDate = endDate,xlim=c(215,(lengthLastYear+220)))
      abline(v=(lengthLastYear-forecastLength+1),col="purple")
      
      ##Add on data points
      points(time.p,p,pch=20,col="red")
      points(time.p,mn,col="green",pch=3)
      points(time.p,me,col="green",pch=1)
      legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
    }

  }
}
dev.off()
