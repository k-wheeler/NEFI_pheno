##' Fills the missing data due to lag in NOAA met with GEFS forecasts
##'
##' @param fileName The GEFS filename
##' @param dataDirectory The directory that the data is stored in
##' @import ncdf4
##' @export
fillNOAAlag <- function(days,siteName){
  dataDirectory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/",siteName,"/",sep="")
  vals <- numeric()
  for(d in 1:length(days)){
    #print(days[d])
    dataDirectoryDate <- paste(dataDirectory,days[d],"/",sep="")
    GEFS_Files <- dir(path=dataDirectoryDate)
    dayVals <- numeric()
    for(f in 1:length(GEFS_Files)){
      #print(GEFS_Files[f])
      weatherFile <- nc_open(paste(dataDirectoryDate,GEFS_Files[f],sep=""))
      Tair <- mean(ncvar_get(weatherFile,"air_temperature")[1:4])-273 ##Convert to C from Kelvin
      nc_close(weatherFile)
      dayVals <- c(dayVals,Tair)
    }
    vals <- rbind(vals,dayVals)
  }
  return(vals)
}
