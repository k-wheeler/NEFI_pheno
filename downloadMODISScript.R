
library("PhenoForecast")
library("PhenologyBayesModeling")
library("coda")
library("dplyr")
library("rjags")
library("MODISTools")
library(doParallel)

##Set and register cores for parallel
n.cores <- 5
registerDoParallel(cores=n.cores)

##Read in data
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 14

endDate <- (Sys.Date()-2)

iseq <- c(seq(1,6),8,9,10,seq(15,27))
output <- 
  foreach(i=iseq) %dopar% {
    #for(i in iseq){
    siteName <- as.character(siteData[i,1])
    print(siteName)

    lat <- as.numeric(siteData[i,2])
    long <- as.numeric(siteData[i,3])
    startDate <- as.Date(siteData[i,7])
    
    mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=startDate,end=endDate,site_name = paste(siteName,"_rel",sep=""),out_dir = dataDirectory,internal=FALSE)
    
    downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
    downloadMODIS(startDate=startDate,endDate=endDate,metric="EVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
    
}
