##' Loads NOAA met data
##'
##' @param station The name of the NOAA station
##' @param startDate
##' @param endDate
##' @import rnoaa
##' @import tidyverse
##' @export
load_NOAA_met <- function(station,startDate="",endDate="") {
  options(noaakey = "fjVseeIJrOLasppGbfwDrYZVsdQaQoCd")
  #station_data <- ghcnd_stations()

  # sites_data <- data.frame(id = siteName,
  #                          latitude = lat,
  #                          longitude = long)
  # nearby_stations <- meteo_nearby_stations(lat_lon_df = sites_data, station_data = station_data, year_min = 2019, year_max = 2019, radius = 50, var = "all")
  # if(siteName=="russellsage"){ #For some reason the closest station to russellsage needs a key code to access the data
  #   station <- as.character(nearby_stations[[1]][2,1])
  # }
  # else if(siteName=="missouriozarks"){
  #   station <- as.character(nearby_stations[[1]][5,1])
  # }
  # else{
  #   station <- as.character(nearby_stations[[1]][1,1])
  # }
  if(nchar(startDate)==0){
    startDate <- as.Date("2019-01-01")
    endDate <- (Sys.Date()-1)
  }
  dat = meteo_tidy_ghcnd(stationid = station,var = c("TAVG","tmin","tmax"), date_min = startDate, date_max = endDate)
  lastDate <- dat$date[length(dat$date)]

  NOAAavgs <- rowMeans(cbind(dat$tmax,dat$tmin))
  NOAAavgs <- NOAAavgs/10 ##Downloads in tenths of a degree C
  for(i in seq((lastDate+1),endDate,"day")){
    NOAAavgs <- c(NOAAavgs,-9999) ## Done to indicate data lag
  }
  return(NOAAavgs)
}
