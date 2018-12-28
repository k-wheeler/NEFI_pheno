library("ncdf4")
library(plyr)

getABI_Index <- function(lat,long,orbitVersion){
  #Function to determine the index of the ABI fixed grid that corresponds to a geodetic (normal) latitude and longitude.
  #lat and long should be in radians not degrees
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

getDataIndex <- function(vals,ch,orbitVersion){
  #Takes in the y,x coordinates on the ABI grid (vals) and returns a list of the i,j index for the product matrices
  #ch refers to the channel number (2,3, or ACM)
  #From PUG_L1b-vol3 pg 15
  #NDVI should use channel 3
  #All add_offset values changed when they moved the satellite
  if(orbitVersion=="OLD"){
    if(ch==3 || ch==5){ ##Ch 3 and Ch 5 have the same scaling numbers (I checked) probably because they are both NIR
      x.scale_factor <- 2.80000003840541e-05 #values from netcdf files under dimensions:x:scale_factor
      x.add_offset <- -0.0750259980559349
      y.scale_factor <- -2.80000003840541e-05
      y.add_offset <- 0.126545995473862
    }
    if(ch==2){
      x.scale_factor <- 1.40000001920271e-05
      x.add_offset <- -0.075033001601696
      y.scale_factor <- -1.40000001920271e-05
      y.add_offset <- 0.126552999019623
    }
    if(ch=="ACM"){
      x.scale_factor <- 5.60000007681083e-05
      x.add_offset <- -0.0750119984149933
      y.scale_factor <- -5.60000007681083e-05
      y.add_offset <- 0.126532003283501
    }
  }
  else{
    if(ch==3 || ch==5){
      x.scale_factor <- 2.80000003840541e-05 #values from netcdf files under dimensions:x:scale_factor
      x.add_offset <- -0.101346001029015
      y.scale_factor <- -2.80000003840541e-05
      y.add_offset <- 0.12822599709034
    }
    if(ch==2){
      x.scale_factor <- 1.40000001920271e-05
      x.add_offset <- -0.101352997124195
      y.scale_factor <- -1.40000001920271e-05
      y.add_offset <- 0.128233000636101
    }
    if(ch=="ACM"){
      x.scale_factor <- 5.60000007681083e-05
      x.add_offset <- -0.101332001388073
      y.scale_factor <- -5.60000007681083e-05
      y.add_offset <- 0.128212004899979
    }
  }
  i <- (vals[2]-x.add_offset)/x.scale_factor
  j <- (vals[1]-y.add_offset)/y.scale_factor
  return(c(i,j))
}

#Test example from PUG_L1b-vol3
#getABI_Index(0.590726966,-1.47813561,orbitVersion="OLD") #passed by returning (0.095340,-0.024052)

#Test for coordinates:
#flor.lat <- 25.204*2*pi/360
#flor.long <- -80.7098*2*pi/360
#getDataIndex(getABI_Index(flor.lat,flor.long),3,orbitVersion="OLD")#NDVI will use channel 3
