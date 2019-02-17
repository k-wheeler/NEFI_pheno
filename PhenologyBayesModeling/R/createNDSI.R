##' Function to create a file with NDSI values for a specific day and time
##'
##' @param day.time The day and time value in the format year,DOY,time (e.g. "20170151200" or "20171821658")
##' @import ncdf4
##' @import plyr
##' @export
createNDSI <- function(day.time){
  #Function to create a file with NDSI values for a specific day/time
  #outputs a csv file
  #day.time needs to be in the format "20171821658"
  #will need to change this to account for the files being in a different folder
  ##NOTE: NOT FINISHED!!!!!!!!!!!
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
