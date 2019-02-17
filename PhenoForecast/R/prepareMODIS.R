#' Function to prepare MODIS data for input into forecast model
#'
#' @param startDate The start date of data in date format
#' @param endDate The end date of data in date format
#' @param metric NDVI or EVI
#' @param timeForecast The time values (in days) based in PhenoCam available data
#' @param dataDirectory The file path for the directory to download and store data in
#' @param lat The latitude for the MODIS in decimal degrees (not necessary if you do not need to download data)
#' @param long The longitude for the MODIS in decimal degrees (not necessary if you do not need to download data)
#' @param siteName The site name for file naming purposes
#' @export
#' @import MODISTools
prepareMODIS <- function(startDate,endDate,metric,timeForecast,dataDirectory,lat="",long="",siteName){
  fileName <- paste(dataDirectory,siteName,"_",metric,"_MOD13Q1_",startDate,"_",endDate,".csv",sep="")
  #print(fileName)
  if(!file.exists(fileName)){ ##If data file is not present, download (Note: This sometimes throws an time out error if it has to be called here)
    downloadMODIS(startDate=startDate,endDate=endDate,metric="NDVI",dataDirectory=dataDirectory,lat=lat,long=long,siteName=siteName)
  }
  ##Read in MODIS data and reformat
  dat <- read.csv(fileName,header=TRUE)
  MODIS.x <- as.Date(dat$calendar_date)
  MODIS.y <- as.numeric(dat$data)/10000

  ##Prepare for phenology forecast
  m <- rep(NA,length(timeForecast))
  for(i in 1:length(MODIS.x)){
    m[which(timeForecast==MODIS.x[i])] <- MODIS.y[i]
  }

  return(m)
}
