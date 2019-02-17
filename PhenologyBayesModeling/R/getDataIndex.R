##' Takes in the y,x coordinates on the ABI grid (vals) and returns a list of the i,j index for the product matrices
##'
##' @param vals The y,x coordinates on the ABI grid
##' @param ch The channel number (2, 3, or ACM). Note: NDVI should use channel 3
##' @export
getDataIndex <- function(vals,ch,orbitVersion){
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
    else if(ch==2){
      x.scale_factor <- 1.40000001920271e-05
      x.add_offset <- -0.075033001601696
      y.scale_factor <- -1.40000001920271e-05
      y.add_offset <- 0.126552999019623
    }
    else if(ch=="ACM"){
      x.scale_factor <- 5.60000007681083e-05
      x.add_offset <- -0.0750119984149933
      y.scale_factor <- -5.60000007681083e-05
      y.add_offset <- 0.126532003283501
    }
  }
  else if(orbitVersion=="NEW"){
    if(ch==3 || ch==5){
      x.scale_factor <- 2.80000003840541e-05 #values from netcdf files under dimensions:x:scale_factor
      x.add_offset <- -0.101346001029015
      y.scale_factor <- -2.80000003840541e-05
      y.add_offset <- 0.12822599709034
    }
    else if(ch==2){
      x.scale_factor <- 1.40000001920271e-05
      x.add_offset <- -0.101352997124195
      y.scale_factor <- -1.40000001920271e-05
      y.add_offset <- 0.128233000636101
    }
    else if(ch=="ACM"){
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


#Test for coordinates:
#flor.lat <- 25.204*2*pi/360
#flor.long <- -80.7098*2*pi/360
#getDataIndex(getABI_Index(flor.lat,flor.long),3,orbitVersion="OLD")#NDVI will use channel 3
