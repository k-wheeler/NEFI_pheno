##' Calculates Normalized Difference Snow Index (NDSI)
##'
##' @param R2 The reflectance from channel 2
##' @param R5 The reflectance from channel 5
calNDSI <- function(R2,R5){
  return((R5-R2)/(R5+R2))
}
