##' Calculates cummulative Tair
##'
##' @param Tairs The start date of desired data download
##' @param days The desired dates
##' @import purrr
##' @export
calSf <- function(Tairs,days) {
  if(Tairs[1]>0){
    Rf <- 28.4/(1+exp(3.4-0.185*Tairs[1]))
  }
  else{
    Rf <- 0
  }
  cumTotal <- Rf
  Sfs <- Rf
  for(i in 2:length(Tairs)){
    if((lubridate::day(days[i])==1)&&(lubridate::month(days[i])==1)){ ##Reset every year
      cumTotal <- 0
    }
    #print(Tairs[i])
    if(!is.na(Tairs[i])){
      if(Tairs[i]>0){
        Rf <- 28.4/(1+exp(3.4-0.185*Tairs[i]))
      }
      else{
        Rf <- 0
      }
    }
    else{
      Rf <- 0
    }
    cumTotal <- cumTotal + Rf
    Sfs <- c(Sfs,cumTotal)
  }
  return(Sfs)
}
