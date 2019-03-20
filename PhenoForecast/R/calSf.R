##' Calculates cummulative Tair
##'
##' @param Tairs The air temperature values in a vector
##' @param dates The desired dates
##' @export
calSf <- function(Tairs,dates) {
  #print(dates)
  Rf <- calRf(Tairs[1])
  cumTotal <- Rf
  Sfs <- Rf
  for(i in 2:length(Tairs)){
    #print(dates[i])
    if((lubridate::day(dates[i])==1)&&(lubridate::month(dates[i])==1)){ ##Reset every year
      cumTotal <- 0
    }
    #print(Tairs[i])
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
