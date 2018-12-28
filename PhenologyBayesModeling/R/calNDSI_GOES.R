#!/usr/bin/env Rscript

library("ncdf4")
library(plyr)

calNDSI <- function(R2,R5){
  return((R5-R2)/(R5+R2))
}

createNDSI <- function(day.time){
  #Function to create a file with NDSI values for a specific day/time
  #outputs a csv file
  #day.time needs to be in the format "20171821658"
  #will need to change this to account for the files being in a different folder
  filestrC02 <- paste("OR_ABI-L1b-RadC-M3C02_G16_s",day.time,sep="")
  filestrC05 <- paste("OR_ABI-L1b-RadC-M3C05_G16_s",day.time,sep="")
  R5 <- ncvar_get(nc_open(dir(pattern=filestrC05)),"Rad")
  R2 <- ncvar_get(nc_open(dir(pattern=filestrC02)),"Rad") #the full R2 dataset
  NDVI.vals <- matrix(ncol=3000,nrow=5000)
  for(i in seq(1,nrow(R2),2)){
    for(j in seq(1,ncol(R2),2)){
      R2.val <- mean(R2[i,j],R2[(i+1),j],R2[i,(j+1)],R2[(i+1),(j+1)])
      R5.val <- R5[(i/2),(j/2)]
      NDSI.vals[(i/2),(j/2)] <- calNDSI(R2.val,R5.val)
    }
    if(i%%500==1){
      print(i) #done to keep track of progress
    }
  }
  file.output.name <- paste("GOES16_NDVI_",day.time,".csv",sep="")
  write.table(NDVI.vals,file.output.name,row.names=FALSE,sep=",")
  return(NDVI.vals)
}

getSpecificNDVI <- function(ind2,ind5,day.time){
  #Function to calculate the NDVI when you have a known pixel
  #print("Entered getSpecificNDVI")
  filestrC02 <- paste("OR_ABI-L1b-RadC-M3C02_G16_s",day.time,sep="")
  filestrC05 <- paste("OR_ABI-L1b-RadC-M3C05_G16_s",day.time,sep="")
  filepathC02 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC02),sep="")
  filepathC05 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC05),sep="")
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
    print("R5 before DQF:")
    print(R5[i5,j5])
    if(!is.na(R5.DQF[i5,j5]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[(i2+1),j2]) && !is.na(R2.DQF[i2,(j2+1)]) && !is.na(R2.DQF[(i2+1),(j2+1)])){
      if(R5.DQF[i5,j5]==0 && R2.DQF[i2,j2]==0 && R2.DQF[i2,j2]==0 && R2.DQF[(i2+1),j2]==0 && R2.DQF[i2,(j2+1)]==0 && R2.DQF[(i2+1),(j2+1)]==0){
        R5.val <- R5[i5,j5]
        R2.val <- mean(R2[i2,j2],R2[(i2+1),j2],R2[i2,(j2+1)],R2[(i2+1),(j2+1)])
        print("")
        print(c(R2.val,R5.val))
        output <- calNDVI(R2.val,R5.val)
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
