##' Calculates Normalized Difference Vegetation Index (NDVI)
##'
##' @param R2 The reflectance from channel 2
##' @param R3 The reflectance from channel 3
calNDVI <- function(R2,R3){
  return((R3-R2)/(R3+R2))
}
