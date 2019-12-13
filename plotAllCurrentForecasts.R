##Plot all autumn forecasts for the current year...
#!/usr/bin/env Rscript
library("PhenoForecast")
library("PhenologyBayesModeling")
library("coda")
library("dplyr")
library("rjags")

##Read in data
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="PhenologyForecastData/"
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 15
dates <- seq(as.Date("2019-08-01"),as.Date("2019-12-10"),"day")
iseq <- c(1,2,3,4,6,15,16,20,18,24)
outFileName <- paste("PhenologyForecastOutput_allSites_allDays_Autumn.pdf",sep="")
pdf(outFileName,height=10,width=6)
par(mfrow=c(5,3))
for(i in iseq){  
  siteName <- as.character(siteData[i,1])
  print(siteName)
  for(d in 1:length(dates)){
    endDate <- dates[d]
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
      lengthLastYear <- (as.numeric(format(endDate,"%j"))+forecastLength)
      print(lengthLastYear)
      ##Random Walk:
      RWFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_",index,"_outBurn.RData",sep="")
      if(file.exists(RWFile)){
        print("Loading randomWalk")
        load(RWFile)
        
        ##Plot Parameter Density Plots
        out.mat.par <- data.frame(as.matrix(outBurnRW$params))
        ##Plot Forecast including calibration
        dayNumber <- dim(as.matrix(outBurnRW$predict))[2]
        out.mat.RW <- as.matrix(outBurnRW$predict)
        
        offset <- 365-as.numeric(format(startDate,"%j"))
        
        ##Plot Current year
        
        lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
        out.mat.lastYear <- out.mat.RW[,lastYearIndices]
        
        plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                           out.mat=out.mat.lastYear,forecastType = "randomWalk",
                           days=seq(1,lengthLastYear,1),endDate=endDate)#,
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
        #legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
      }
      ##Basic logistic:
      lengthLastYear <- lengthLastYear - 212 ##Done for autumn
      LFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_basicLog_",index,"_outBurn.RData",sep="")
      if(file.exists(LFile)){
        print("Loading basicLog")
        load(LFile)
        out.mat.par <- data.frame(as.matrix(outBurnL$params))
        print(colnames(out.mat.par))
        
        ##Plot
        dayNumber <- dim(as.matrix(outBurnL$predict))[2]
        out.mat.L <- as.matrix(outBurnL$predict)
        
        ##Plot Current year
        
        lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
        out.mat.lastYear <- out.mat.L[,lastYearIndices]
        
        plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,
                           out.mat=out.mat.lastYear,forecastType = "Logistic",days=seq(213,(lengthLastYear+212),1),
                           endDate=endDate,xlim=c(215,(lengthLastYear+220)))#,
        
        abline(v=(lengthLastYear-forecastLength+1),col="purple")
        
        ##Add on data points
        points(time.p,p,pch=20,col="red")
        points(time.p,mn,col="green",pch=3)
        points(time.p,me,col="green",pch=1)
        #legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
      }
      ##Logistic with covariates model
      LCFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_CDD_20_meanTair_outBurn_GCC.RData",sep="")
      print(LCFile)
      if(file.exists(LCFile)){
        print("Loading CDD20")
        load(LCFile)
        outBurnLC <- outBurn
        out.mat.par <- data.frame(as.matrix(outBurnLC$params))
        ##Plot
        dayNumber <- dim(as.matrix(outBurnLC$predict))[2]
        out.mat.LC <- as.matrix(outBurnLC$predict)
        
        ##Plot Current year
        lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
        out.mat.lastYear <- out.mat.LC[,lastYearIndices]
        
        plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.lastYear,
                           forecastType = "CDD_20",days=seq(213,(lengthLastYear+212),1),
                           endDate = endDate,xlim=c(215,(lengthLastYear+220)))
        abline(v=(lengthLastYear+212-forecastLength+1),col="purple")
        
        ##Add on data points
        points(time.p,p,pch=20,col="red")
        points(time.p,mn,col="green",pch=3)
        points(time.p,me,col="green",pch=1)
        #legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))
      }else{
        plot(seq(1,lengthLastYear,1),rep(NA,lengthLastYear),ylim=c(0,1))
      }
    }
  }
}
dev.off()
