###Script to determine met stations
library("PhenoForecast")
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="PhenologyForecastData/"
endDate <- Sys.Date()-1
startDate <- as.Date("2019-07-01")
iseq <- c(1,2,3,4,6,15,16,20,18,24)
#iseq <- c(2,3,4,6,15,16,20,18,24)
s <- 1
stations <- character()
sites <- character()
for(s in iseq){
  siteName <- as.character(siteData[s,1])
  print(siteName)
  
  lat <- as.numeric(siteData[s,2])
  long <- as.numeric(siteData[s,3])
  #startDate <- as.Date(siteData[s,7])
  elev <- as.numeric(siteData$elev[s])
  
  station <- identifyNOAAstation(lat = lat,long = long,elev = elev,
                                 siteName = siteName,startDate = startDate,endDate = endDate)
  print(station)
  stations <- c(stations,station)
  sites <- c(sites,siteName)
}
output <- cbind(sites,stations)
write.table(output,file="PhenologyForecastData/met_stations.csv",col.names = FALSE,row.names = FALSE,sep = ",")
