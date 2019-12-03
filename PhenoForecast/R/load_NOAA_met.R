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

  if(nchar(startDate)==0){
    startDate <- as.Date("2019-08-01")
    endDate <- (Sys.Date()-1)
  }
  dat <- getStationData(station=station, startDate=startDate,endDate=endDate)
  #dat <-  meteo_tidy_ghcnd(stationid = station,var = c("TAVG","tmin","tmax"), date_min = startDate, date_max = endDate)

  dat[is.na(dat$tmax),]$tmin <- -99990
  dat[is.na(dat$tmax),]$tmax <- -99990
  lastDate <- dat$date[length(dat$date)]

  NOAAavgs <- rowMeans(cbind(dat$tmax,dat$tmin))
  NOAAavgs <- NOAAavgs/10 ##Downloads in tenths of a degree C

  if(lastDate<endDate){
    missingDays <- seq((lastDate+1),endDate,"day")
    #print(missingDays)
    if(length(missingDays)>0){
      for(i in 1:length(missingDays)){
        NOAAavgs <- c(NOAAavgs,-9999) ## Done to indicate data lag
      }
    }
  }
  return(NOAAavgs)
}
