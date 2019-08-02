###Forecast simulation for Harvard Forest

library("PhenoForecast")
library("PhenologyBayesModeling")
library("coda")
library("dplyr")
library("rjags")
##Need to make individual forecast plots
#setwd("/Users/Kathryn/Documents/PhD_Research/NEFI_pheno")
#siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
#dataDirectory="PhenologyForecastData/"
i <- 1 ##For Harvard
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 15
mdlStr <- "LC3"

siteName <- "harvard"
forecastFiles <- intersect(dir(path=paste(dataDirectory,"ForecastOutputs/AllForecasts",sep=""),pattern=siteName),
                           dir(path=paste(dataDirectory,"ForecastOutputs/AllForecasts",sep=""),pattern=mdlStr))
for(f in 1:length(forecastFiles)){
  endDate <- as.Date(strsplit(forecastFiles[f],"_")[[1]][3])
  print(endDate)
  lengthLastYear <- (as.numeric(format(endDate,"%j"))+forecastLength)
  siteName <- as.character(siteData[i,1])
  
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
  phenoData <- phenoData[phenoData$date<endDate,]
  
  p <- phenoData$gcc_mean[phenoData$year==2019]
  
  time.p.date <- as.Date(phenoData$date[phenoData$year==2019])
  time.p <-  as.numeric(format(time.p.date,"%j"))
  
  p <- rescale(c=cMeans.p[length(cMeans.p)],d=dMeans.p[length(dMeans.p)],yseq=p)
  p[p<0] <- 0
  
  ##MODIS NDVI and EVI
  mn <- prepareMODIS(startDate=startDate,endDate=endDate,metric="NDVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
  me <- prepareMODIS(startDate=startDate,endDate=endDate,metric="EVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
  mn <- rescale(yseq=mn,c=cMeans.mn[length(cMeans.mn)],d=dMeans.mn[length(dMeans.mn)])
  me <- rescale(yseq=me,c=cMeans.me[length(cMeans.me)],d=dMeans.me[length(dMeans.me)])
  mn[mn<0] <- 0
  me[me<0] <- 0
  load(paste(dataDirectory,"ForecastOutputs/AllForecasts/",forecastFiles[f],sep=""))
  
  #outBurnLC <- outBurnLC2
  out.mat.par <- data.frame(as.matrix(outBurnLC$params))
  colnames(out.mat.par)
  
  ##Plot
  par(mfrow=c(1,1))
  dayNumber <- dim(as.matrix(outBurnLC$predict))[2]
  out.mat.LC <- as.matrix(outBurnLC$predict)
  
  ##Plot Current year
  lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
  out.mat.lastYear <- out.mat.LC[,lastYearIndices]
  plotFileName <- paste("ForecastPlots/",siteName,"_",endDate,"_",mdlStr,".png",sep="")
  png(file=plotFileName, width=10, height=5,units="in",res=1000)
  plotForecastOutput(siteName=siteName,URL=URL[length(URL)],forecastLength=forecastLength,
                     out.mat=out.mat.lastYear,forecastType = "Logistic Covariate",
                     days=seq(1,lengthLastYear,1),xlim=c(10,160),
                     plotTitle = "Harvard Forest Forecast",endDate=endDate)
  abline(v=(lengthLastYear-forecastLength+1),col="purple")
  
  ##Add on data points
  points(time.p,p,pch=20,col="red")
  points(time.p,mn,col="blue",pch=3,cex=2)
  points(time.p,me,col="blue",pch=1,cex=2)
  legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","blue","blue"),pch=c(20,3,1))
  dev.off()

}






