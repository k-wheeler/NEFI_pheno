##' Calculates cummulative Tair
##'
##' @param Tairs The air temperature values in a vector
##' @param dates The desired dates
##' @export
calSf <- function(Tairs,dates) {
  Rf <- calRf(Tairs[1])
  cumTotal <- Rf
  Sfs <- Rf
  for(i in 2:length(dates)){
    if((lubridate::day(dates[i])==1)&&(lubridate::month(dates[i])==1)){ ##Reset every year
      cumTotal <- 0
    }
    if(!is.na(Tairs[i])){
      Rf <- calRf(Tairs[i])
    }
    else{
      Rf <- 0
    }
    cumTotal <- cumTotal + Rf
    Sfs <- c(Sfs,cumTotal)
  }
  return(Sfs)
}
