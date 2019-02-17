##' Loads the GEFS weather forecast data
##'
##' @param fileName The GEFS filename
##' @import ncdf4
##' @export
load_GEFS_Forecast<- function(fileName) {
  weatherFile <- nc_open(fileName)
  Tair <- ncvar_get(weatherFile,"air_temperature")
  Tairs <- numeric()
  for(i in seq(1,length(Tair),4)){
    Tairs <- c(Tairs,(mean(Tair[i:(i+3)])-273))
  }
  return(Tairs)
}
