##' Calculates the index of the ABI fixed grid that corresponds to a geodetic (normal) latitute and longitude
##'
##' @param lat Latitude in radians (not degrees)
##' @param long Longitude in radians (not degrees)
##' @param orbitVersion The orbit version ("OLD" if before the satellite moved into its final orbit for GOES-East)
##' @import ncdf4
##' @import plyr
##' @export
getABI_Index <- function(lat,long,orbitVersion="NEW"){
  #Note: to index into the ABI fixed grid, you index by [y,x] not [x,y].
  #From PUG_L1b-vol3 pg 23

  #values from netcdf metadata:
  H <- 42164160  #goes_imagery_projection:perspective_point_height + goes_imagery_proj:semi-major_axis
  if(orbitVersion=="OLD"){
    long0 <- -1.562069680 #goes_imagery_projection:longitude_of_projection_origin
  }
  else{
    long0 <- -1.308996939 #goes_imagery_projection:longitude_of_projection_origin
  }
  r.pol <- 6356752.31414     #goes_imagery_projection:semi_minor_axis
  r.eq <- 6378137 #goes_imagery_projection:semi_major_axis
  e.value <- 0.0818191910435
  lat.c <- atan(r.pol**2/r.eq**2*tan(lat)) #geocentric latitude
  rc <- r.pol/(sqrt(1-e.value**2*(cos(lat.c))**2)) #geocentric distance to the point on the ellipsoid
  Sx <- H-rc*cos(lat.c)*cos(long-long0)
  Sy <- -rc*cos(lat.c)*sin(long-long0)
  Sz <- rc*sin(lat.c)
  y <- atan(Sz/Sx) #N/S Elevation Angle
  x <- asin(-Sy/(sqrt(Sx**2+Sy**2+Sz**2))) #E/W Scanning Angle
  return(c(y,x))
}

#Test example from PUG_L1b-vol3
#getABI_Index(0.590726966,-1.47813561,orbitVersion="OLD") #passed by returning (0.095340,-0.024052)
