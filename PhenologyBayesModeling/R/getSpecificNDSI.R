##' Function to calculate the NDSI when you have a known pixel
##'
##' @param day.time The day and time value in the format year,DOY,time (e.g. "20170151200" or "20171821658")
##' @param ind2 The index for the site in the channel 2 grid
##' @param ind5 The index for the site in the channel 5 grid
##' @param dataFolder The data folder where the data is located (e.g. "GOES_Data2017/")
##' @import ncdf4
##' @import plyr
##' @export
getSpecificNDSI <- function(ind2,ind5,day.time,dataFolder="GOES_Data2017/"){
  #Function to calculate the NDSI when you have a known pixel
  #print("Entered getSpecificNDVI")
  filestrC02 <- paste("OR_ABI-L1b-RadC-M3C02_G16_s",day.time,sep="")
  filestrC05 <- paste("OR_ABI-L1b-RadC-M3C05_G16_s",day.time,sep="")
  filepathC02 <- paste(dataFolder,dir(path=dataFolder,pattern=filestrC02),sep="")
  filepathC05 <- paste(dataFolder,dir(path=dataFolder,pattern=filestrC05),sep="")
  if(file.exists(filepathC02) && file.exists(filepathC05) && !dir.exists(filepathC02) && !dir.exists(filepathC05)){
    R2.file <-nc_open(filepathC02)
    R5.file <-nc_open(filepathC05)

    R5 <- ncvar_get(R5.file,"Rad")
    R2 <- ncvar_get(R2.file,"Rad") #the full R2 dataset
    R5.kappa0 <- ncvar_get(R5.file,"kappa0")
    R2.kappa0 <- ncvar_get(R2.file,"kappa0")
    R5.DQF <- ncvar_get(R5.file,"DQF") #Data Quality Flags
    R2.DQF <- ncvar_get(R2.file,"DQF")
    R5 <- R5 * R5.kappa0 #done to covert radiance to reflectance
    R2 <- R2 * R2.kappa0

    i2 <- ind2[1]
    j2 <- ind2[2]
    i5 <- ind5[1]
    j5 <- ind5[2]
    if(!is.na(R5.DQF[i5,j5]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[(i2+1),j2]) && !is.na(R2.DQF[i2,(j2+1)]) && !is.na(R2.DQF[(i2+1),(j2+1)])){
      if(R5.DQF[i5,j5]==0 && R2.DQF[i2,j2]==0 && R2.DQF[i2,j2]==0 && R2.DQF[(i2+1),j2]==0 && R2.DQF[i2,(j2+1)]==0 && R2.DQF[(i2+1),(j2+1)]==0){
        R5.val <- R5[i5,j5]
        R2.val <- mean(R2[i2,j2],R2[(i2+1),j2],R2[i2,(j2+1)],R2[(i2+1),(j2+1)])
        output <- calNDSI(R2.val,R5.val)
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
  return(output)
}
