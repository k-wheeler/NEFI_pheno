##' Loads the GEFS weather forecast data
##'
##' @param fileName The GEFS filename
##' @param dataDirectory The directory that the data is stored in
##' @import ncdf4
##' @export
load_GEFS_Forecast<- function(fileName,dataDirectory) {
  ##Extract info on the start date and time
  startDateTime <- strsplit(fileName,split="[.]")[[1]][4]
  startTime <- strsplit(startDateTime,"T")[[1]][2]
  startDate <-strsplit(startDateTime,"T")[[1]][1]

  ##Open the file
  weatherFile <- nc_open(paste(dataDirectory,fileName,sep=""))
  Tair <- ncvar_get(weatherFile,"air_temperature")
  timeVals <- weatherFile$dim$time$vals
  dates <- as.POSIXlt(paste(startDate,startTime,sep=" "))+(timeVals*3600) ##Convert to date and times
  currentDate <- lubridate::date(dates[1])

  ##Average out the Tairs that have the same UTC day
  Tairs <- numeric()
  subTair <- numeric()
  for(i in 1:length(dates)){
    if(lubridate::date(dates[i])==currentDate){
      subTair <- c(subTair,Tair[i])
    }
    else{
      Tairs <- c(Tairs,mean(subTair))
      currentDate <- lubridate::date(dates[i])
      subTair <- Tair[i]
    }

  }
  Tairs <- c(Tairs,mean(subTair))
  if(startTime!="00:00"){ ##Removes the partial start and end dates
    Tairs <- Tairs[2:(length(Tairs)-1)]
  }
  Tairs <- Tairs-273 ##Convert to C from Kelvin
  return(Tairs)
}
