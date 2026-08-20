##' Calculates Rf (from Melaas et al. 2016 GCB)
##'
##' @param Tair The air temperature
##' @export
calRf <- function(Tair){
  if(!is.na(Tair)){
    if(Tair>0){
      Rf <- 28.4/(1+exp(3.4-0.185*Tair))
    }
    else{
      Rf <- 0
    }
  }
  else{
    Rf <- 0
  }
}
