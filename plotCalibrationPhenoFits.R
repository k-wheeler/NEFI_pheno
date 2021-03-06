#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)
library(PhenoForecast)
library(PhenologyBayesModeling)
library(coda)
library(rjags)
library(doParallel)
library(ecoforecastR)
library("MODISTools")

#season <- "spring"
#endDate <- (Sys.Date()-1)
#startDate <- as.Date("2013-01-01")
endDate <- as.Date("2020-01-27")
forecastLength <- 0
iseq <- c(2,4,5,seq(21,27))
iseq <- seq(1,27)
iseq <- 17
iseq <- c(1,2,3,4,6,15,16,19,20,22)
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

for(i in iseq){
  if(as.character(siteData$PFT[i])=="DB"){
    siteName <- as.character(siteData[i,1])
    print(siteName)

    lat <- as.numeric(siteData[i,2])
    long <- as.numeric(siteData[i,3])
    startDate <- (as.Date(siteData[i,7]))
    URL <- as.character(siteData$URL[i])
    URL2 <- as.character(siteData$URL2[i])
    URL3 <- as.character(siteData$URL3[i])
    if(nchar(URL2)>0){
      URL <- c(URL,URL2)
      if(nchar(URL3)>0){
        URL <- c(URL,URL3)
      }
    }
    days <- seq(as.Date(startDate),(as.Date(endDate)+forecastLength),"day")
    dataDirectory="PhenologyForecastData/"
    
    ##Download/load data
    ##Download new MODIS data
    #downloadMODIS(startDate=startDate,endDate=endDate,metric="rel",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
    
    #downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
    #downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
    
    ##PhenoCam data
    newMonths <- lubridate::month(days)
    newYears <- lubridate::year(days)
    
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
    years <- seq(lubridate::year(startDate),2019)
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
    
     pdf(file=paste(siteName,"_PhenologyForecast_previousFitsNEW.pdf",sep=""),height=6,width=10)
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
        print(colnames(var.mat))
        cMeans.p <- c(cMeans.p,mean(var.mat.DF$c))
        dMeans.p <- c(dMeans.p,mean(var.mat.DF$d))
        ks <- var.mat.DF$k
        if(length(ks)==0){
          ks <- c(182,182)
        }
        kMeans.p <- c(kMeans.p,mean(ks))
        
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
        print(colnames(var.mat))
        cMeans.mn <- c(cMeans.mn,mean(var.mat.DF$c))
        dMeans.mn <- c(dMeans.mn,mean(var.mat.DF$d))
        ks <- var.mat.DF$k
        if(length(ks)==0){
          ks <- c(182,182)
        }
        kMeans.mn <- c(kMeans.mn,mean(ks))
        
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
        print(colnames(var.mat))
        cMeans.me <- c(cMeans.me,mean(var.mat.DF$c))
        dMeans.me <- c(dMeans.me,mean(var.mat.DF$d))
        ks <- var.mat.DF$k
        if(length(ks)==0){
          ks <- c(182,182)
        }
        kMeans.me <- c(kMeans.me,mean(ks))
        
        CI <- createCI(PFT="DB",var.mat = var.mat,xseq=DOYs,doRescale = FALSE,seasonOrder = "SF")
        plot(DOYs,me.yr,pch=20,main=paste(years[j],"ME"))
        ciEnvelope(x=DOYs,ylo=CI[1,],yhi=CI[3,],col="lightblue")
        points(DOYs,me.yr,pch=20)
      }
    dev.off()
    ##Write files of c, d, and k means
    write.table(cbind(cMeans.p,dMeans.p,kMeans.p,years),row.names = FALSE,col.names = TRUE,file=paste("PhenologyForecastData/",siteName,"_forecast_phenoFits_PC2020.csv",sep=""),sep=",")
    write.table(cbind(cMeans.mn,dMeans.mn,kMeans.mn,years),row.names = FALSE,col.names = TRUE,file=paste("PhenologyForecastData/",siteName,"_forecast_phenoFits_MN2020.csv",sep=""),sep=",")
    write.table(cbind(cMeans.me,dMeans.me,kMeans.me,years),row.names = FALSE,col.names = TRUE,file=paste("PhenologyForecastData/",siteName,"_forecast_phenoFits_ME2020.csv",sep=""),sep=",")
  }
}
