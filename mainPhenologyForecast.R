install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)
library("PhenoForecast")
library("PhenologyBayesModeling")
library("coda")
library("dplyr")
library("rjags")
#library(doParallel)

##Set and register cores for parallel
#n.cores <- 6
#registerDoParallel(cores=n.cores)

##Read in data
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="PhenologyForecastData/"
forecastLength <- 16

#endDate <- (Sys.Date()-1)
#startDate <- as.Date("2013-01-01")
endDate <- as.Date("2019-01-27")
i <- 10
#Create Forecast outputs
#output <- 
#  foreach(i=1:nrow(siteData)) %dopar% {
siteName <- as.character(siteData[i,1])
print(siteName)
GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/willowcreek/2019-01-25/",sep="")
GEFS_files <- dir(path=GEFS_Directory,pattern="NOAA_GEFS")

print(GEFS_files)
URL <- as.character(siteData[i,4])
lat <- as.numeric(siteData[i,2])
long <- as.numeric(siteData[i,3])
startDate <- as.Date(siteData[i,7])
##Download new MODIS data
downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)

##Load rescaling data
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

##Create Random Walk forecast if needed
# outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
# if(!file.exists(outputFile)){
#   outBurnRW <- phenologyForecast(forecastType = "randomWalk",forecastLength = forecastLength,siteName=siteName,URL=URL,lat=lat,long=long,dataDirectory=dataDirectory,startDate,endDate,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me)
#   if(typeof(outBurn)!=typeof(FALSE)){
#     save(outBurnRW,file=outputFile)
#   }
# }
# 
# ##Create logistic forecast if needed
# outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_logistic_outBurn.RData",sep="")
# if(!file.exists(outputFile)){
#   outBurnL <- phenologyForecast(forecastType = "logistic",forecastLength = forecastLength,siteName=siteName,URL=URL,lat=lat,long=long,dataDirectory=dataDirectory,startDate=startDate,endDate=endDate,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me)
#   if(typeof(outBurn)!=typeof(FALSE)){
#     save(outBurn,file=outputFile)
#   }
# }
##Create a Logistic with Covariate Model
if(siteName=="willowCreek"){
  outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_LC_outBurn.RData",sep="")
  if(!file.exists(outputFile)){
    #########################
    forecastLength=16
    outBurnLC <- phenologyForecast(forecastType = "logisticCov",forecastLength = forecastLength,siteName="willowCreek",URL=URL,lat=lat,long=long,dataDirectory=dataDirectory,as.Date(startDate),as.Date(endDate),GEFS_Files=GEFS_files,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me,GEFS_Directory = GEFS_Directory)
    save(outBurnLC,file=outputFile)
    #   variableNames <- c("p.PC","p.MN","p.ME","p.proc","x","b0","b1")
    #   
    #   
    #   
    #   out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames)
    #   plot(dat$p,dat$Sf,pch=20)
    #   if(typeof(outBurn)!=typeof(FALSE)){
    #     save(outBurn,file=outputFile)
    #   }
  }
}

# rescaleObs()
# 
# 
# ##Download meteological data
# Tair <- download_US_WCr_met(start_date=as.Date("2016-01-01"), end_date=as.Date("2018-06-30"))
# 
# days <- seq(startDate,endDate,"day")
# newTairs <- rep(NA,length(days))
# for(i in 1:nrow(dat)){
#   newTairs[which(days==dat[i,1])] <- dat[i,4]
# }
# newMonths <- lubridate::month(days)
# newYears <- lubridate::year(days)
# dat2 <- data.frame(dates=days,years=newYears,months=newMonths,Tairs=newTairs)
# 
# ##Need to separate out spring data (months January-June)
# URL <- "https://phenocam.sr.unh.edu/data/archive/willowcreek/ROI/willowcreek_DB_1000_1day.csv"
# phenoData <- download.phenocam(URL)
# p.old <- phenoData$gcc_mean
# time.old <-  as.Date(phenoData$date)
# p <- rep(NA,length(days))
# for(i in 1:length(p.old)){
#   p[which(days==time.old[i])] <- p.old[i]
# }
# 
# 
# dat2 <- dat2[dat2$months%in%seq(1,6,1),]
# 
# startDate <- as.Date("2016-01-01")
# endDate <- as.Date("2018-06-30")
# startDate <- as.Date("2012-04-26")
# 
# Tair <- matrix(nrow=181,ncol=0)
# for(yr in lubridate::year(startDate):lubridate::year(endDate)){
#   subDat <- dat2[dat2$years==yr,]
#   print(length(subDat$Tairs))
#   Tair <- cbind(Tair,subDat$Tairs)
# }
# 
# dates <- as.Date(character())
# current.date <- as.Date(dat$date[1])
# Tairs <- numeric()
# sub.Tairs <- numeric()
# years <- numeric()
# months <- numeric()
# for(i in 1:nrow(dat)){
#   print(i)
#   if(as.Date(dat$date[i])==as.Date(current.date)){
#     if(dat$Tair[i]!=-999){
#       sub.Tairs <- c(sub.Tairs,dat$Tair[i])
#     }
#   }
#   else{
#     Tairs <- c(Tairs,mean(sub.Tairs))
#     dates <- c(dates,current.date)
#     months <- c(months,dat$Month)
#     years <- c(years,dat[(i-1)]$Year)
#     current.date <- as.Date(dat$date[i])
#     sub.Tairs <- numeric()
#   }
#   
# }
# 
# }
#plot(dates,Tairs,pch=20)

##Plotting Forecasts
# pdf(file=paste("ForecastOutputs_",endDate,".pdf",sep=""),width=12,height=8)
# for(i in 1:16){
#   siteName <- as.character(siteData[i,1])
#   print(siteName)
#   URL <- as.character(siteData[i,4])
#   lat <- as.numeric(siteData[i,2])
#   long <- as.numeric(siteData[i,3])
#   startDate <- as.Date(siteData[i,7])
#   outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
#   load(outputFile)
# days <- c(dat2$dates,rep(NA,forecastLength))
# dayNumber <- dim(as.matrix(outBurnRW$predict))[2]
# plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,outBurn=outBurnRW,forecastType = "randomWalk",days=seq(1,dayNumber,1))#2484
# points(days,rep(1,length(days)),pch=20)
# offset <- 365-as.numeric(format(startDate,"%j"))
# for(i in seq(offset,2484,365)){
#   abline(v=i,col="red")
# }
# 
# dayNumber <- dim(as.matrix(outBurnL$predict))[2]
# plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,outBurn=outBurnL,forecastType = "Logistic",days=seq(1,dayNumber,1))#2484
# points(days,rep(1,length(days)),pch=20)
# #offset <- 365-as.numeric(format(startDate,"%j"))
# for(i in seq(1,2484,182)){
#   abline(v=i,col="red")
# }

#   
#   outputFile <- paste("PhenologyForecastData/",siteName,"_",startDate,"_",endDate,"_logistic_outBurn.RData",sep="")
#   load(outputFile)
#   plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,dataSource="PhenoCam",outBurn=outBurn,forecastType = "logistic")
# }
# dev.off()
