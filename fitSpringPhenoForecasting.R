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
registerDoParallel(cores=n.cores)

iseq <- c(seq(2,6),8,9,11,seq(15,20))
#iseq <- c(10)
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

for(i in 1:nrow(siteData)){
  if(as.character(siteData$PFT[i])=="DB"){
    siteName <- as.character(siteData[i,1])
    print(siteName)
    URL <- as.character(siteData[i,4])
    lat <- as.numeric(siteData[i,2])
    long <- as.numeric(siteData[i,3])
    startDate <- (as.Date(siteData[i,7]))
    URL2 <- as.character(siteData$URL2[i])
    URL3 <- as.character(siteData$URL3[i])
    startDate2 <- as.character(siteData$startDate2[i])
    days <- seq(as.Date(startDate),(as.Date(endDate)),"day")
    dataDirectory="PhenologyForecastData/"
    
    ##Download/load data
    ##Download new MODIS data
    lastDate <- (as.Date(startDate) - 1)
    newDQFFileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new DQF data downloaded
    if(!file.exists(newDQFFileName)){
      print("Downloading MODIS DQF File")
      try(mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=(lastDate+1),end=endDate,site_name = paste(siteName,"_rel",sep=""),out_dir = dataDirectory,internal=FALSE),silent=TRUE)
      downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
      downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
    }
    
    ##PhenoCam data
    newMonths <- lubridate::month(days)
    newYears <- lubridate::year(days)
    #print("Done with newYears")

    phenoData <- download.phenocam(URL) 
    if(nchar(URL2>0)){
      phenoData2 <- download.phenocam(URL2)
      phenoData3 <- rbind(phenoData,phenoData2)
      phenoData3 <- phenoData3[order(phenoData3[,'date']),]
      phenoData <- phenoData3[!duplicated(phenoData3$date),]
    }
    if(nchar(URL3>0)){
      phenoData2 <- download.phenocam(URL3)
      phenoData3 <- rbind(phenoData,phenoData2)
      phenoData3 <- phenoData3[order(phenoData3[,'date']),]
      phenoData <- phenoData3[!duplicated(phenoData3$date),]
    }
    
    p.old <- phenoData$gcc_mean
    time.old <-  as.Date(phenoData$date)
    p <- rep(NA,length(days))
    for(i in 1:length(p.old)){
      p[which(days==time.old[i])] <- p.old[i]
    }
    dat2 <- data.frame(dates=days,years=newYears,months=newMonths,p=p)
    dat2$mn <- prepareMODIS(startDate=startDate,endDate=endDate,metric="NDVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)
    dat2$me <- prepareMODIS(startDate=startDate,endDate=endDate,metric="EVI",timeForecast=days,dataDirectory=dataDirectory,siteName=siteName)
    
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
    
    output <- 
      foreach(j=1:length(years)) %dopar% {
        # for(j in 1:length(years)){
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
        ks <- var.mat.DF$k
        if(length(ks)==0){
          ks <- c(182,182)
        }
        kMeans.p <- c(kMeans.p,mean(ks))
        
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
        ks <- var.mat.DF$k
        if(length(ks)==0){
          ks <- c(182,182)
        }
        kMeans.mn <- c(kMeans.mn,mean(ks))
        
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
        ks <- var.mat.DF$k
        if(length(ks)==0){
          ks <- c(182,182)
        }
        kMeans.me <- c(kMeans.me,mean(ks))

      }
 }
}
