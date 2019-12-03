#' Calculates a met station score based on Hopkin's Bioclimatic Law
#'
#' @param siteLat
#' @param siteLong
#' @param siteElev
#' @param statLat
#' @param statLong
#' @param statElev
#'
#' @return
#' @export
#'
#' @examples
calStationScore <- function(siteLat,siteLong, siteElev, statLat,statLong, statElev){
  ##Based on Hopkin's Bioclimatic Law
  latDif <- abs(siteLat-statLat) #Degrees
  longDif <- abs(siteLong-statLong) #Degrees
  elevDif <- abs(siteElev-statElev) * 3.2808 #Convert meeters to feet
  overallScore <- 100 - latDif - longDif/5 - elevDif/400
  return(overallScore)
}
