#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)
library(PhenoForecast)
library(PhenologyBayesModeling)
library(coda)
library(rjags)
library(doParallel)
library(ecoforecastR)

endDate <- as.Date("2020-01-01")
forecastLength <- 0

##Set and register cores for parallel
n.cores <- 8
registerDoParallel(cores=n.cores)

siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

foreach(i=1:nrow(siteData)) %dopar% {
#for(i in 1:nrow(siteData)){
  siteName <- as.character(siteData[i,1])
  print(siteName)
  
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
  startDate <- (as.Date(siteData[i,7]))
  #startDate <- as.Date("2017-01-01")
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

  for(j in 1:length(p.old)){
    p[which(days==time.old[j])] <- p.old[j]

  }
  dat2 <- data.frame(dates=days,years=newYears,months=newMonths,p=p)
  data2 <- data.frame(p=dat2$p)
  
  p <- matrix(nrow=365,ncol=0)
  
  years <- seq(lubridate::year(startDate),2019)
  for(j in years){
    subDat <- data2[lubridate::year(as.Date(dat2$dates))==j,]
    p <- cbind(p,subDat)
    
  }
  #dataFinal <- list(p=p,mn=mn,me=me)
  DOYs <- seq(1,365,1)
  
  DB.vars <- c("TranF","bF","TranS","bS","c","d","prec")
  
  cMeans.p <- numeric()
  dMeans.p <- numeric()
  kMeans.p <- numeric()
  
  pdf(file=paste(siteName,"_PhenologyForecast_previousFits_GCC.pdf",sep=""),height=6,width=10)
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
    
  }
  dev.off()
  ##Write files of c, d, and k means
  write.table(cbind(cMeans.p,dMeans.p,kMeans.p,years),row.names = FALSE,col.names = TRUE,file=paste("PhenologyForecastData/",siteName,"_forecast_phenoFits_PC2020.csv",sep=""),sep=",")
}

