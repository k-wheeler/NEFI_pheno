install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)
library(PhenoForecast)
library(PhenologyBayesModeling)
library(coda)
library(rjags)
library(doParallel)
##Create Phenology Fits for willow Creek spring data (should actually fit spring and autumn together)
#season <- "spring"
#endDate <- (Sys.Date()-1)
#startDate <- as.Date("2013-01-01")
endDate <- as.Date("2019-01-27")
forecastLength <- 0

n.cores <- 6

#register the cores.
#registerDoParallel(cores=n.cores)

i <- 10
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
siteName <- as.character(siteData[i,1])
print(siteName)
URL <- as.character(siteData[i,4])
lat <- as.numeric(siteData[i,2])
long <- as.numeric(siteData[i,3])
startDate <- as.Date(siteData[i,7])
days <- seq(as.Date(startDate),(as.Date(endDate)+forecastLength),"day")
dataDirectory="PhenologyForecastData/"


##Download/load data
##Download new MODIS data
#downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
#downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)

##PhenoCam data
newMonths <- lubridate::month(days)
newYears <- lubridate::year(days)
#print("Done with newYears")
PC.fileName <- paste(dataDirectory,siteName,"_",startDate,"_",endDate,"_PC_Data.RData",sep="")
if(!file.exists(PC.fileName)){
  phenoData <- download.phenocam(URL) #<- doesn't work on the cluster I think (so need to download all of the phenodata and save it)
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
for(i in (lubridate::year(as.Date(dat2$dates[1]))+1):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){
  subDat <- data2[lubridate::year(as.Date(dat2$dates))==i,]
  p <- cbind(p,subDat$p)
  mn <- cbind(mn,subDat$mn)
  me <- cbind(me,subDat$me)
}
#dataFinal <- list(p=p,mn=mn,me=me)
DOYs <- seq(1,365,1)
##Divide the data up into each spring
#DB.vars <- c("TranS","bS","c","d","prec","k")
DB.vars <- c("TranF","bF","TranS","bS","c","d","prec","k")
#cMeans <- numeric()
#dMeans <- numeric()
#kMeans <- numeric()
j=1
years <- seq(2013,2018)
#output <- 
#  foreach(j=1:6) %dopar% {
for(j in 4:6){
  ##PhenoCam Fits
  outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_PC_",years[j],"_varBurn.RData",sep="")
  if(!file.exists(outFileName)){
  p.yr <- p[,j]
  #plot(DOYs,p.yr,pch=20)
  data <- list(x=DOYs,y=p.yr,n=length(p.yr))
  j.model <- createModel_DB(data=data,dataSource = "PC.GCC",seasonOrder = "SF")
  varBurn <- runMCMC_Model(j.model = j.model,variableNames = DB.vars,baseNum=40000,iterSize=20000)
  #out.mat <- as.matrix(varBurn)
  #cMeans <- c(cMeans,mean(out.mat[,5]))
  #dMeans <- c(dMeans,mean(out.mat[,6]))
  #kMeans <- c(kMeans,mean(out.mat[,7]))

  save(varBurn,file=outFileName)
  }
  
  ##MODIS NDVI Fits
  outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_MN_",years[j],"_varBurn.RData",sep="")
  if(!file.exists(outFileName)){
    mn.yr <- mn[,j]
    data <- list(x=DOYs,y=mn.yr,n=length(mn.yr))
    j.model <- createModel_DB(data=data,dataSource = "MODIS.NDVI",seasonOrder = "SF")
    varBurn <- runMCMC_Model(j.model = j.model,variableNames = DB.vars,baseNum=40000,iterSize=20000)
    save(varBurn,file=outFileName)
  }
  
  ##MODIS EVI Fits
  outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_ME_",years[j],"_varBurn.RData",sep="")
  if(!file.exists(outFileName)){
  me.yr <- me[,j]
  data <- list(x=DOYs,y=me.yr,n=length(me.yr))
  j.model <- createModel_DB(data=data,dataSource = "MODIS.EVI",seasonOrder = "SF")
  varBurn <- runMCMC_Model(j.model = j.model,variableNames = DB.vars,baseNum=40000,iterSize=20000)
  save(varBurn,file=outFileName)
  }
  
}
#plot(DOYs,p.yr,pch=20)
#years <- c(2013,2014,2015,2016,2017,2018)
#write.table(cbind(cMeans,dMeans,kMeans,years),row.names = FALSE,col.names = TRUE,file="willowCreek_forecast_phenoFits.csv",sep=",")

##Run phenology fits