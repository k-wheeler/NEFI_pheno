#' Function to rescale phenology measurements to percent canopy cover for that calendar year
#'
#' @param times A vector of days (x-values)
#' @param vals A vector of the measurements corresponding for the days (y-values)
#' @param cVals The c values for rescaling
#' @param dVals The d values for rescaling
#' @param partialStart Whether or not there is a partial year of data included (e.g. with random walk)
#' @import PhenologyBayesModeling
#' @export
rescaleObs <- function(times,vals,cVals,dVals,partialStart=FALSE,yr=2019){
  new.ys <- numeric()
  years <- seq((lubridate::year(times[1])),yr)
  if(partialStart){
    yr <- 1
    year.vals <- vals[(lubridate::year(times)==(years[yr]))] ##The observations for first partial year
    year.c <- cVals[1]
    year.d <- dVals[1]
    year.vals.new <- rescale(c=year.c,d=year.d,yseq=year.vals)
    new.ys <- c(new.ys,year.vals.new)
    for(yr in 2:(length(years)-1)){
      #print(years[yr])
      year.vals <- vals[(lubridate::year(times)==years[yr])] ##The observations for each year
      if(years[yr]==2018){
        #print(year.vals)
      }
      year.c <- cVals[(yr-1)]
      year.d <- dVals[(yr-1)]

      year.vals.new <- rescale(c=year.c,d=year.d,yseq=year.vals)
      new.ys <- c(new.ys,year.vals.new)
      #print(range(year.vals,na.rm = TRUE))
      #print(range(year.vals.new,na.rm=TRUE))
    }
  }
  else{
    for(yr in 1:(length(years)-1)){
      year.vals <- vals[(lubridate::year(times)==years[yr])] ##The observations for each year
      year.c <- cVals[yr]
      year.d <- dVals[yr]
      year.vals.new <- rescale(c=year.c,d=year.d,yseq=year.vals)
      new.ys <- c(new.ys,year.vals.new)
    }
  }

  yr <- length(years)
  year.vals <- vals[(lubridate::year(times)==years[yr])] ##The observations for the new year
  year.c <- cVals[length(cVals)]
  year.d <- dVals[length(dVals)]
  year.vals.new <- rescale(c=year.c,d=year.d,yseq=year.vals)
  new.ys <- c(new.ys,year.vals.new)
  #plot(times,new.ys,pch=20,main="Inside rescale")
  # new.ys <- numeric()
  # for(yr in ((lubridate::year(times[1])+1):2018)){
  #   year.vals <- vals[(lubridate::year(times)==yr)] ##The observations for each year
  #   year.max <- max(year.vals,na.rm=TRUE) ##Max for that year
  #   year.min <- min(year.vals,na.rm=TRUE)  ##Min for that year
  #   year.vals.new <- (year.vals-year.max)/(year.max-year.min)+1 ##Rescale so the range is between (0,1)
  #   new.ys <- c(new.ys,year.vals.new)
  # }
  # ##Rescale partial start year
  # yr <- lubridate::year(range(times)[1])
  # year.vals <- vals[(lubridate::year(times)==yr)]
  # next.year.vals <- vals[(lubridate::year(times)==(yr+1))]
  # year.max <- max(next.year.vals,na.rm=TRUE)
  # year.min <- min(next.year.vals,na.rm=TRUE)
  # year.vals.new <- (year.vals-year.max)/(year.max-year.min)+1
  # new.ys <- c(year.vals.new,new.ys)
  #
  # ##Rescale partial current year
  # yr <- 2019
  # year.vals <- vals[(lubridate::year(times)==yr)]
  # past.year.vals <- vals[(lubridate::year(times)==(yr-1))]
  # year.max <- max(past.year.vals,na.rm=TRUE)
  # year.min <- min(past.year.vals,na.rm=TRUE)
  # year.vals.new <- (year.vals-year.max)/(year.max-year.min)+1
  # new.ys <- c(new.ys,year.vals.new)
  return(new.ys)
}
