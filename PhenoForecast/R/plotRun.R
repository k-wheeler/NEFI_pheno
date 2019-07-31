#' plotRun
#'
#' @param out.mat
#' @param forecastType
#' @param endDate
#' @return
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#' @export
#'
plotRun <- function(out.mat,forecastType,endDate){
  ##Plot Current year
  dayNumber <- ncol(out.mat)
  #endDate <- as.Date("2019-01-23")
  lengthLastYear <- (as.numeric(format(endDate,"%j")))
  #print(lengthLastYear)
  lastYearIndices <- seq(((dayNumber-lengthLastYear)+1),dayNumber,1)
  out.mat.lastYear <- out.mat[,lastYearIndices]

  plotForecastOutput(siteName=siteName,URL=URL,forecastLength=0,out.mat=out.mat.lastYear,forecastType = forecastType,days=seq(1,lengthLastYear,1),xlim=c(0,175),endDate=endDate)


}
