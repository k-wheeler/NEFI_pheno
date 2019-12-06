##' Identifies the best NOAA station
##'
##' @param lat Site latitude
##' @param long Site longitude
##' @param elev Site elevation (m)
##' @param startDate Start date
##' @param endDate End Date
##' @import rnoaa
##' @import tidyverse
##' @export
identifyNOAAstation <- function(lat,long,elev,siteName,startDate="", endDate=""){
  options(noaakey = "fjVseeIJrOLasppGbfwDrYZVsdQaQoCd")
  station_data <- ghcnd_stations()

  sites_data <- data.frame(id = siteName,
                           latitude = lat,
                           longitude = long)
  nearby_stations <- meteo_nearby_stations(lat_lon_df = sites_data, station_data = station_data,
                                           year_min = 2019, year_max = 2019, radius = 100, var = "all")
  if(nchar(startDate)==0){
    startDate <- as.Date("2019-08-01")
    endDate <- (Sys.Date()-1)
  }
  stationScores <- numeric()
  scoredStations <- character()
  for(i in 1:nrow(nearby_stations[[1]])){
    print(i)
    station <- as.character(nearby_stations[[1]][i,1])
    dat <- getStationData(station=station, startDate=startDate,endDate=endDate)
    if(typeof(dat)!=typeof(FALSE)){
      if(nrow(dat)>0){
        if(sum(is.na(dat$tmax))<5){
          statLat <- as.numeric(nearby_stations[[1]][i,3])
          statLong <- as.numeric(nearby_stations[[1]][i,4])
          statElev <- as.numeric(station_data[which(station_data$id==station)[1],]$elevation)
          stationScores <- c(stationScores,calStationScore(siteLat=lat,siteLong=long,
                                                           siteElev = elev ,
                                                           statLat=statLat,
                                                           statLong = statLong ,
                                                           statElev = statElev))
          scoredStations <- c(scoredStations,station)
        }
      }
    }
  }
  station <- scoredStations[which.max(stationScores)]
  return(station)
}




