#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)
library(PhenoForecast)
library(PhenologyBayesModeling)
library(coda)
library(rjags)
library(doParallel)
library(ecoforecastR)
library("MODISTools")
library("doParallel")
##Create Phenology Fits for willow Creek spring data (should actually fit spring and autumn together)
#season <- "spring"
#endDate <- (Sys.Date()-1)
#startDate <- as.Date("2013-01-01")
endDate <- as.Date("2019-01-27")
forecastLength <- 0

n.cores <- 6

#register the cores.
#registerDoParallel(cores=n.cores)

#iseq <- c(seq(1,6),8,9,11,seq(15,20))
iseq <- c(1)
for(i in iseq){
  siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
  siteName <- as.character(siteData[i,1])
  print(siteName)
  URL <- as.character(siteData[i,4])
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
  startDate <- (as.Date(siteData[i,7]))
  days <- seq(as.Date(startDate),(as.Date(endDate)+forecastLength),"day")
  dataDirectory="PhenologyForecastData/"
  
  
  ##Download/load data
  ##Download new MODIS data
  lastDate <- (as.Date(startDate) - 1)
  newDQFFileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new DQF data downloaded
  if(!file.exists(newDQFFileName)){
    print("Downloading MODIS DQF File")
    try(mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=(lastDate+1),end=endDate,site_name = paste(siteName,"_rel",sep=""),out_dir = dataDirectory,internal=FALSE),silent=TRUE)
  }

  downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
  downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)

  ##PhenoCam data
  newMonths <- lubridate::month(days)
  newYears <- lubridate::year(days)
  #print("Done with newYears")
  PC.fileName <- paste(dataDirectory,siteName,"_",startDate,"_",endDate,"_PC_Data.RData",sep="")
  if(!file.exists(PC.fileName)){
    phenoData <- download.phenocam(URL) 
    save(phenoData,file=PC.fileName)
  }
  load(PC.fileName)
  
  p.old <- phenoData$gcc_mean
  time.old <-  as.Date(phenoData$date)
  p <- rep(NA,length(days))
  for(i in 1:length(p.old)){
    p[which(days==time.old[i])] <- p.old[i]
  }
  dat2 <- data.frame(dates=days,years=newYears,months=newMonths,p=p)
  dat2$mn <- prepareMODIS(startDate=startDate,endDate=endDate,metric="NDVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)
  dat2$me <- prepareMODIS(startDate=startDate,endDate=endDate,metric="EVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)
  #print("Done with MODIS")
  #print(dim(dat2))
  #dat2 <- dat2[dat2$months%in%seq(1,7,1),]
  #print(dim(dat2))
  
  data2 <- data.frame(p=dat2$p)
  
  data2$mn <- dat2$mn
  data2$me <- dat2$me
  
  p <- matrix(nrow=365,ncol=0)
  mn <- matrix(nrow=365,ncol=0)
  me <- matrix(nrow=365,ncol=0)
  years <- seq(lubridate::year(startDate),2018)
  for(i in years){
    subDat <- data2[lubridate::year(as.Date(dat2$dates))==i,]
    p <- cbind(p,subDat$p)
    mn <- cbind(mn,subDat$mn)
    me <- cbind(me,subDat$me)
  }
  #dataFinal <- list(p=p,mn=mn,me=me)
  DOYs <- seq(1,365,1)
  ##Divide the data up into each spring
  #DB.vars <- c("TranS","bS","c","d","prec","k")
  DB.vars <- c("TranF","bF","TranS","bS","c","d","prec")
  #cMeans <- numeric()
  #dMeans <- numeric()
  #kMeans <- numeric()
  j=1
  #years <- seq(2013,2018)
  cMeans.p <- numeric()
  dMeans.p <- numeric()
  kMeans.p <- numeric()
  cMeans.mn <- numeric()
  dMeans.mn <- numeric()
  kMeans.mn <- numeric()
  cMeans.me <- numeric()
  dMeans.me <- numeric()
  kMeans.me <- numeric()
  
  #pdf(file=paste(siteName,"PhenologyForecast_previousFitsNEW.pdf",sep=""),height=6,width=10)
  output <- 
    foreach(j=1:length(years)) %dopar% {
      for(j in 1:length(years)){
      print(years[j])
      ##PhenoCam Fits
      outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_PC_",years[j],"_varBurn.RData",sep="")
      p.yr <- p[,j]
      if(!file.exists(outFileName)){
        data <- list(x=DOYs,y=p.yr,n=length(p.yr))
        j.model <- createModel_DB(data=data,dataSource = "PC.GCC",seasonOrder = "SF")
        varBurn <- runMCMC_Model(j.model = j.model,variableNames = DB.vars,baseNum=40000,iterSize=20000)
        save(varBurn,file=outFileName)
      }
      load(outFileName)
      var.mat <- as.matrix(varBurn)
      var.mat.DF <- data.frame(var.mat)
      cMeans.p <- c(cMeans.p,mean(var.mat.DF$c))
      dMeans.p <- c(dMeans.p,mean(var.mat.DF$d))
      kMeans.p <- c(kMeans.p,mean(var.mat.DF$k))
      
      CI <- createCI(PFT="DB",var.mat = var.mat,xseq=DOYs,doRescale = FALSE,seasonOrder = "SF")
      plot(DOYs,p.yr,pch=20,main=paste(years[j],"PC"))
      ciEnvelope(x=DOYs,ylo=CI[1,],yhi=CI[3,],col="lightblue")
      points(DOYs,p.yr,pch=20)
      
      
      ##MODIS NDVI Fits
      outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_MN_",years[j],"_varBurn.RData",sep="")
      mn.yr <- mn[,j]
      if(!file.exists(outFileName)){
        data <- list(x=DOYs,y=mn.yr,n=length(mn.yr))
        j.model <- createModel_DB(data=data,dataSource = "MODIS.NDVI",seasonOrder = "SF")
        varBurn <- runMCMC_Model(j.model = j.model,variableNames = DB.vars,baseNum=40000,iterSize=20000)
        save(varBurn,file=outFileName)
      }
      load(outFileName)
      var.mat <- as.matrix(varBurn)
      var.mat.DF <- data.frame(var.mat)
      cMeans.mn <- c(cMeans.mn,mean(var.mat.DF$c))
      dMeans.mn <- c(dMeans.mn,mean(var.mat.DF$d))
      kMeans.mn <- c(kMeans.mn,mean(var.mat.DF$k))
      
      CI <- createCI(PFT="DB",var.mat = var.mat,xseq=DOYs,doRescale = FALSE,seasonOrder = "SF")
      plot(DOYs,mn.yr,pch=20,main=paste(years[j],"MN"))
      ciEnvelope(x=DOYs,ylo=CI[1,],yhi=CI[3,],col="lightblue")
      points(DOYs,mn.yr,pch=20)
      
      ##MODIS EVI Fits
      outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_ME_",years[j],"_varBurn.RData",sep="")
      me.yr <- me[,j]
      if(!file.exists(outFileName)){
        data <- list(x=DOYs,y=me.yr,n=length(me.yr))
        j.model <- createModel_DB(data=data,dataSource = "MODIS.EVI",seasonOrder = "SF")
        varBurn <- runMCMC_Model(j.model = j.model,variableNames = DB.vars,baseNum=40000,iterSize=20000)
        save(varBurn,file=outFileName)
      }
      load(outFileName)
      var.mat <- as.matrix(varBurn)
      var.mat.DF <- data.frame(var.mat)
      cMeans.me <- c(cMeans.me,mean(var.mat.DF$c))
      dMeans.me <- c(dMeans.me,mean(var.mat.DF$d))
      kMeans.me <- c(kMeans.me,mean(var.mat.DF$k))
      
      CI <- createCI(PFT="DB",var.mat = var.mat,xseq=DOYs,doRescale = FALSE,seasonOrder = "SF")
      plot(DOYs,me.yr,pch=20,main=paste(years[j],"ME"))
      ciEnvelope(x=DOYs,ylo=CI[1,],yhi=CI[3,],col="lightblue")
      points(DOYs,me.yr,pch=20)
    }
  #dev.off()
  
  ##Write files of c, d, and k means
  write.table(cbind(cMeans.p,dMeans.p,kMeans.p,years),row.names = FALSE,col.names = TRUE,file=paste("PhenologyForecastData/",siteName,"_forecast_phenoFits_PC.csv",sep=""),sep=",")
  write.table(cbind(cMeans.mn,dMeans.mn,kMeans.mn,years),row.names = FALSE,col.names = TRUE,file=paste("PhenologyForecastData/",siteName,"_forecast_phenoFits_MN.csv",sep=""),sep=",")
  write.table(cbind(cMeans.me,dMeans.me,kMeans.me,years),row.names = FALSE,col.names = TRUE,file=paste("PhenologyForecastData/",siteName,"_forecast_phenoFits_ME.csv",sep=""),sep=",")
  
}
