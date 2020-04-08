#'Calculate NDVI for a specific day and time
#'
#' @param siteData matrix of site data where the sites are in individual rows and the columns are siteName, latitude, longitude, timezone
#' @param day.time The year,day,time (e.g. 20172001020 for 10:20 on the 200th day of 2017)
#' @param orbitVersion NEW or OLD depending on the position of the satellite (only days before 17 November 2017 are OLD). Default value is NEW
#' @param dataPath The directory where all of the GOES data is located
#' @param mVersion The mode/version of the ABI (3,4, or 6)
#' @export
#' @import ncdf4
#' @import plyr
createNDVI_sub <- function(siteData,day.time,orbitVersion=NEW,dataPath,mVersion){
  print("Entered createNDVI_sub")
  ##Create File Paths
  ACM.path <- paste(dataPath,"/",dir(path=dataPath,pattern=paste("OR_ABI-L2-ACMC-M",mVersion,"_G16_s",day.time,sep="")),sep="")
  filestrC03 <- paste("OR_ABI-L1b-RadC-M",mVersion,"C03_G16_s",day.time,sep="")
  filestrC02 <- paste("OR_ABI-L1b-RadC-M",mVersion,"C02_G16_s",day.time,sep="")
  filePathC02 <- paste(dataPath,"/",dir(path=dataPath,pattern=filestrC02),sep="")
  filePathC03 <- paste(dataPath,"/",dir(path=dataPath,pattern=filestrC03),sep="")
  NDVIs <- numeric()

  ##Open Files
  if(nchar(filePathC02)>20 && nchar(filePathC03)>20){
    if(file.exists(filePathC02) && file.exists(filePathC03) && !dir.exists(filePathC02) && !dir.exists(filePathC03)){
      ACM.file <-nc_open(ACM.path)
      R2.file <- nc_open(filePathC02)
      R3.file <- nc_open(filePathC03)

      ##Extract Data
      R3 <- ncvar_get(R3.file,"Rad")
      R2 <- ncvar_get(R2.file,"Rad") #the full R2 dataset
      R3.kappa0 <- ncvar_get(R3.file,"kappa0")
      R2.kappa0 <- ncvar_get(R2.file,"kappa0")
      R3.DQF <- ncvar_get(R3.file,"DQF") #Data Quality Flags
      R2.DQF <- ncvar_get(R2.file,"DQF")
      R3 <- R3 * R3.kappa0 #done to covert radiance to reflectance
      R2 <- R2 * R2.kappa0
      clouds <- ncvar_get(ACM.file,"BCM")
      clouds.DQF <- ncvar_get(ACM.file,"DQF")

      for(i in 1:nrow(siteData)){
        ##General Site Data
        siteName <- as.character(siteData[i,1])
        lat <- as.numeric(siteData[i,2])
        long <- as.numeric(siteData[i,3])
        TZ <- as.numeric(siteData[i,4])

        ##Determine index values
        lat.rd <- lat*2*pi/360
        long.rd <- long*2*pi/360

        Ind2 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion=orbitVersion),2,orbitVersion=orbitVersion)
        Ind3 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion=orbitVersion),3,orbitVersion=orbitVersion)
        ACM.ind <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion=orbitVersion),"ACM",orbitVersion=orbitVersion)

        i2 <- Ind2[1]
        j2 <- Ind2[2]
        i3 <- Ind3[1]
        j3 <- Ind3[2]

        if(!is.na(R3.DQF[i3,j3]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[i2,j2]) &&
           !is.na(R2.DQF[(i2+1),j2]) && !is.na(R2.DQF[i2,(j2+1)]) && !is.na(R2.DQF[(i2+1),(j2+1)]) &&
           !is.na(clouds[ACM.ind[1],ACM.ind[2]]) && !is.na(clouds.DQF[ACM.ind[1],ACM.ind[2]])){

          if(R3.DQF[i3,j3]==0 && R2.DQF[i2,j2]==0 && R2.DQF[i2,j2]==0 && R2.DQF[(i2+1),j2]==0 &&
             R2.DQF[i2,(j2+1)]==0 && R2.DQF[(i2+1),(j2+1)]==0 && clouds[ACM.ind[1],ACM.ind[2]] == 0 &&
             clouds.DQF[ACM.ind[1],ACM.ind[2]] == 0){

            R3.val <- R3[i3,j3]
            R2.val <- mean(R2[i2,j2],R2[(i2+1),j2],R2[i2,(j2+1)],R2[(i2+1),(j2+1)])
            output <- calNDVI(R2.val,R3.val)
          }
          else{
            output <- NA
          }
        }
        else{
          output <- NA
        }
        print(paste("NDVI:",output))
        NDVIs <- c(NDVIs,output)
        nc_close(ACM.file)
        nc_close(R2.file)
        nc_close(R3.file)
      }
    }
    else{
      return(rep(NA,nrow(siteData)))
    }
  }
  else{
    return(rep(NA,nrow(siteData)))
  }
  return(NDVIs)
}
