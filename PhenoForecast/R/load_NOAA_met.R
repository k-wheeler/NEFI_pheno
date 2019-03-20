##' Calculates cummulative Tair (within one year for one ensemble)
##'
##' @param lat The site latitude
##' @param long The site longitude
##' @param years The desired years
##' @import rnoaa
##' @import tidyverse
##' @export
load_NOAA_met <- function(lat,long,years) {
  options(noaakey = "fjVseeIJrOLasppGbfwDrYZVsdQaQoCd")
  station_data <- ghcnd_stations()

  sites_data <- data.frame(id = siteName,
                           latitude = lat,
                           longitude = long)
  nearby_stations <- meteo_nearby_stations(lat_lon_df = sites_data, station_data = station_data, year_min = 2018, year_max = 2018, radius = 30, var = "all")
  station <- as.character(nearby_stations[[1]][1,1])
  data = meteo_tidy_ghcnd(stationid = station,var = c("TAVG","tmin","tmax"), date_min = "2019-01-01", date_max = (Sys.Date()-1))

  NOAAavgs <- rowMeans(cbind(dat$tmax,dat$tmin))
  return(NOAAavgs)
}
