#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)
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

endDate <- (Sys.Date()-2)
lengthLastYear <- (as.numeric(format(endDate,"%j"))+forecastLength)
print(lengthLastYear)

iseq <- c(seq(1,6),8,9,10,seq(15,27))
iseq <- c(seq(1,4),6,8,9,10,seq(15,27))
iseq <- c(seq(1,4),6,8,9,10,15,16,seq(18,27))
outFileName <- paste("PhenologyForecastOutput_allSites_LC2_",endDate,".pdf",sep="")
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
    
    ##Load Data
    ##PhenoCam
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
    print(mn)
    print(me)
    
    ##Logistic with covariates model
    LCFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_LC2_outBurn.RData",sep="")
    if(file.exists(LCFile)){
      load(LCFile)
      outBurnLC <- outBurnLC2
      out.mat.par <- data.frame(as.matrix(outBurnLC$params))
      colnames(out.mat.par)
      
      ##Plot
      par(mfrow=c(1,1))
      dayNumber <- dim(as.matrix(outBurnLC$predict))[2]
      out.mat.LC <- as.matrix(outBurnLC$predict)
      
      ##Plot Current year
      lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
      out.mat.lastYear <- out.mat.LC[,lastYearIndices]
      
      plotForecastOutput(siteName=siteName,URL=URL[length(URL)],forecastLength=forecastLength,out.mat=out.mat.lastYear,forecastType = "Logistic Covariate",days=seq(1,lengthLastYear,1))
      abline(v=(lengthLastYear-forecastLength+1),col="purple")
      
      ##Add on data points
      points(time.p,p,pch=20,col="red")
      points(time.p,mn,col="blue",pch=3)
      points(time.p,me,col="blue",pch=1)
      legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","blue","blue"),pch=c(20,3,1))
    }
  }
}


dev.off()
