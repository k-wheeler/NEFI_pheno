##' Calculates cummulative Tair (within one year for one ensemble)
##'
##' @param lat The site latitude
##' @param long The site longitude
##' @param years The desired years
##' @param siteName The name of the site
##' @import rnoaa
##' @import tidyverse
##' @export
load_NOAA_met <- function(lat,long,years,siteName) {
  options(noaakey = "fjVseeIJrOLasppGbfwDrYZVsdQaQoCd")
  options(key="fjVseeIJrOLasppGbfwDrYZVsdQaQoCd")
  station_data <- ghcnd_stations()

  sites_data <- data.frame(id = siteName,
                           latitude = lat,
                           longitude = long)
  nearby_stations <- meteo_nearby_stations(lat_lon_df = sites_data, station_data = station_data, year_min = 2018, year_max = 2018, radius = 50, var = "all")
  if(siteName=="russellsage" || siteName=="missouriozarks"){ #For some reason the closest station to russellsage needs a key code to access the data
    station <- as.character(nearby_stations[[1]][2,1])
  }
  else{
    station <- as.character(nearby_stations[[1]][1,1])
  }
  dat = meteo_tidy_ghcnd(stationid = station,var = c("TAVG","tmin","tmax"), date_min = "2019-01-01", date_max = (Sys.Date()-1))

  NOAAavgs <- rowMeans(cbind(dat$tmax,dat$tmin))
  NOAAavgs <- NOAAavgs/10 ##Downloads in tenths of a degree C
  return(NOAAavgs)
}
