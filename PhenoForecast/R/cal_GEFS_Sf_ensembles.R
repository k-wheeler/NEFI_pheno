##' Loads the GEFS weather forecast data ensembles
##'
##' @param fileNames A list of the GEFS filename
##' @param dates The list of dates for the forecast
##' @import ncdf4
##' @export
cal_GEFS_Sf_ensembles<- function(fileNames) {
  ##Load the Tairs for each ensemble
  for(j in 1:length(fileNames)){
  weatherFile <- nc_open(fileName)
  Tair <- ncvar_get(weatherFile,"air_temperature")
  Tairs <- numeric()
  for(i in seq(1,length(Tair),4)){
    Tairs <- c(Tairs,(mean(Tair[i:(i+3)])-273))
  }
  ##Calculate Sf

  }
  #return(Tairs)
}
