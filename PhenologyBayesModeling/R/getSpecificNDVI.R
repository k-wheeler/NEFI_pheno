##' Function to calculate the NDVI when you have a known pixel
##'
##' @param day.time The day and time value in the format year,DOY,time (e.g. "20170151200" or "20171821658")
##' @param ind2 The index for the site in the channel 2 grid
##' @param ind5 The index for the site in the channel 5 grid
##' @param dataFolder The data folder where the data is located (e.g. "GOES_Data2017/")
##' @import ncdf4
##' @import plyr
##' @export
getSpecificNDVI <- function(ind2,ind3,day.time){
  #Function to calculate the NDVI when you have a known pixel
  #print("Entered getSpecificNDVI")
  filestrC02 <- paste("OR_ABI-L1b-RadC-M3C02_G16_s",day.time,sep="")
  filestrC03 <- paste("OR_ABI-L1b-RadC-M3C03_G16_s",day.time,sep="")
  filepathC02 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC02),sep="")
  filepathC03 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC03),sep="")
  if(file.exists(filepathC02) && file.exists(filepathC03) && !dir.exists(filepathC02) && !dir.exists(filepathC03)){
    R2.file <-nc_open(filepathC02)
    R3.file <-nc_open(filepathC03)

    R3 <- ncvar_get(R3.file,"Rad")
    R2 <- ncvar_get(R2.file,"Rad") #the full R2 dataset
    R3.kappa0 <- ncvar_get(R3.file,"kappa0")
    R2.kappa0 <- ncvar_get(R2.file,"kappa0")
    R3.DQF <- ncvar_get(R3.file,"DQF") #Data Quality Flags
    R2.DQF <- ncvar_get(R2.file,"DQF")
    R3 <- R3 * R3.kappa0 #done to covert radiance to reflectance
    R2 <- R2 * R2.kappa0

    i2 <- ind2[1]
    j2 <- ind2[2]
    i3 <- ind3[1]
    j3 <- ind3[2]
    print("R3 before DQF:")
    print(R3[i3,j3])
    if(!is.na(R3.DQF[i3,j3]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[(i2+1),j2]) && !is.na(R2.DQF[i2,(j2+1)]) && !is.na(R2.DQF[(i2+1),(j2+1)])){
      if(R3.DQF[i3,j3]==0 && R2.DQF[i2,j2]==0 && R2.DQF[i2,j2]==0 && R2.DQF[(i2+1),j2]==0 && R2.DQF[i2,(j2+1)]==0 && R2.DQF[(i2+1),(j2+1)]==0){
        R3.val <- R3[i3,j3]
        R2.val <- mean(R2[i2,j2],R2[(i2+1),j2],R2[i2,(j2+1)],R2[(i2+1),(j2+1)])
        print("")
        print(c(R2.val,R3.val))
        output <- calNDVI(R2.val,R3.val)
      }
      else{
        output <- NA
      }
    }
    else{
      output <- NA
    }
  }
  else{
    output <- NA
  }
  print("Inside function NDVI:")
  print(output)
  return(output)
}
