#' Gets data for a specific station with a tryCatch
#'
#' @param station
#' @param startDate
#' @param endDate
#' @import rnoaa
#' @import tidyverse
#'
#' @return
#' @export
#'
#' @examples
getStationData <- function(station,startDate,endDate){
  #options(noaakey = "fjVseeIJrOLasppGbfwDrYZVsdQaQoCd")
  options(noaakey = "RjtyiLudyeOuOpHOyRnGJnfKVFTFGrTP")
  out <- tryCatch(
    { meteo_tidy_ghcnd(stationid = station,var = c("TAVG","tmin","tmax"),
                       date_min = startDate, date_max = endDate)
    },
    error=function(cond){
      message(station)
      message("Needs a key")
      return(FALSE)
    },
    warning=function(cond){
      message("Here's the original warning message:")
      message(cond)
      return(NULL)
    }
  )
}
