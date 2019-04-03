##Install latest versions of packages
#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)

library("PhenoForecast")
library("PhenologyBayesModeling")
library("coda")
library("dplyr")
library("rjags")
library("MODISTools")
#library(doParallel)

##Set and register cores for parallel
#n.cores <- 6
#registerDoParallel(cores=n.cores)

##Read in data
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

# if (!dir.exists(dataDirectory)) {
#   dir.create(dataDirectory, recursive=TRUE, showWarnings = FALSE)
# }

forecastLength <- 14

endDate <- (Sys.Date()-1)
iseq <- c(1,10)
#Create Forecast outputs
#output <- 
#  foreach(i=1:nrow(siteData)) %dopar% {
for(i in iseq){
  siteName <- as.character(siteData[i,1])
  print(siteName)
  #GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/willowcreek/",endDate,"/",sep="")
  GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/harvardbarn/",endDate,"/",sep="")
  GEFS_files <- dir(path=GEFS_Directory,pattern="NOAA_GEFS")
  
  URL <- as.character(siteData[i,4])
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
  startDate <- as.Date(siteData[i,7])
  ##Download new MODIS data
  ##Download DQF file if there are no previous ones 
  files <- intersect(dir(path=dataDirectory,pattern=paste(siteName,"_rel",sep="")),dir(path=dataDirectory,pattern="MOD13Q1")) #Current downloaded data files
  
  if(length(files)==0){
    newDQFFileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",startDate,"_",endDate,".csv",sep="") #File name for new DQF data downloaded
    if(!file.exists(newDQFFileName)){
      print("Downloading MODIS DQF File because no files are present")
      mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=startDate,end=endDate,site_name = paste(siteName,"_rel",sep=""),out_dir = dataDirectory,internal=FALSE)
    }
  }
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
  
  saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
  
  ##Create Random Walk forecast if needed
  outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
  if(!file.exists(outputFile)){
    outBurnRW <- phenologyForecast(forecastType = "randomWalk",forecastLength = forecastLength,siteName=siteName,URL=URL,lat=lat,long=long,dataDirectory=dataDirectory,startDate,endDate,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me)
    if(typeof(outBurnRW)!=typeof(FALSE)){
      save(outBurnRW,file=outputFile)
    }
  }
  
  ##Create logistic forecast if needed
  outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_logistic_outBurn.RData",sep="")
  if(!file.exists(outputFile)){
    outBurnL <- phenologyForecast(forecastType = "logistic",forecastLength = forecastLength,siteName=siteName,URL=URL,lat=lat,long=long,dataDirectory=dataDirectory,startDate=startDate,endDate=endDate,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me)
    if(typeof(outBurnL)!=typeof(FALSE)){
      save(outBurnL,file=outputFile)
    }
  }
  
  ##Create a Logistic with Covariate Model
  outputFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_LC_outBurn.RData",sep="")
  if(!file.exists(outputFile)){
    outBurnLC <- phenologyForecast(forecastType = "logisticCov",forecastLength = forecastLength,siteName=siteName,URL=URL,lat=lat,long=long,dataDirectory=dataDirectory,as.Date(startDate),as.Date(endDate),GEFS_Files=GEFS_files,cValsPC=cMeans.p,dValsPC=dMeans.p,cValsMN=cMeans.mn,dValsMN=dMeans.mn,cValsME=cMeans.me,dValsME=dMeans.me,GEFS_Directory = GEFS_Directory)
    if(typeof(outBurnLC)!=typeof(FALSE)){
      save(outBurnLC,file=outputFile)
    }
  }
}

