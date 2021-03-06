
#!/usr/bin/env Rscript

library("ncdf4")
library(plyr)
library("PhenologyBayesModeling")
library(doParallel)

#detect cores.
#n.cores <- detectCores()
n.cores <- 5

#register the cores.
registerDoParallel(cores=n.cores)


createNDVI_GOES_diurnal <- function(lat,long,siteID,startDay,endDay){
  #load/calcuate GOES NDVI data
  lat.rd <- lat*2*pi/360
  long.rd <- long*2*pi/360
  
  Ind2 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion="OLD"),2,orbitVersion="OLD")
  Ind3 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion="OLD"),3,orbitVersion="OLD")
  ACM.ind <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion="OLD"),"ACM",orbitVersion="OLD")
  
  i2 <- Ind2[1]
  j2 <- Ind2[2]
  i3 <- Ind3[1]
  j3 <- Ind3[2]
  
  NDVI.vals <- list()
  
  days <- seq(startDay,endDay,1)
  day.time.vals <- list()
  for(i in 1:length(days)){
    #print(days)
    days[i] <- as.numeric(days[i])
    if(as.numeric(days[i]) < 10){
      days[i] <- paste("00",as.character(days[i]),sep="")
    }
    else if(as.numeric(days[i]) < 100){
      days[i] <- paste("0",as.character(days[i]),sep="")
    }
  }
  for (i in 1:length(days)){
    #print(days)
    #days[i] <- as.numeric(days[i])
    #if(days[i] < 10){
    #  days[i] <- paste("00",as.character(days[i]),sep="")
    #}
    #else if(days[i] < 100){
    #  days[i] <- paste("0",as.character(days[i]),sep="")
    #}
    #print(days[i])
    days[i] <- as.character(days[i])
    filestrACM <- paste("OR_ABI-L2-ACMC-M3_G16_s2017",days[i],sep="")
    ACM.files <- dir(path="GOES_Data2017",pattern=filestrACM)
    print(length(ACM.files))
    for(j in 1:length(ACM.files)){
      day.time <- substr(ACM.files[j],24,34)
      #print(j)
      print(day.time)
      day.time.vals <- c(day.time.vals,day.time)
      filePath <- paste("GOES_Data2017/",ACM.files[j],sep="")
      #print(filePath)
      ACM.file <-nc_open(paste("GOES_Data2017/",ACM.files[j],sep=""))
      #print(dim(ncvar_get(ACM.file, "BCM")))
      #print(ACM.ind)
      clouds <- ncvar_get(ACM.file,"BCM")[ACM.ind[1],ACM.ind[2]]
      if(!is.na(clouds)){
        if (clouds ==0){
          filestrC03 <- paste("OR_ABI-L1b-RadC-M3C03_G16_s",day.time,sep="")
          filestrC02 <- paste("OR_ABI-L1b-RadC-M3C02_G16_s",day.time,sep="")
          filePathC02 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC02),sep="")
          filePathC03 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC03),sep="")
          if(nchar(filePathC02)>20 & nchar(filePathC03)>20){
            R2.file <- nc_open(paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC02),sep=""))
            R3.file <- nc_open(paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC03),sep=""))
            R3.DQF <- ncvar_get(R3.file,"DQF")
            R2.DQF <- ncvar_get(R2.file,"DQF")
            if(R3.DQF[i3,j3]==0 & R2.DQF[i2,j2]==0 & R2.DQF[i2,j2]==0 & R2.DQF[(i2+1),j2]==0 & R2.DQF[i2,(j2+1)]==0 & R2.DQF[(i2+1),(j2+1)]==0){
              NDVI.val <- getSpecificNDVI(Ind2,Ind3,day.time)
            }
            else{
	          NDVI.val <- NA
            }
          }
          else{
            NDVI.val <- NA
          }
        }
        else{
          NDVI.val <- NA
        }
        }
      else{
        NDVI.val <- NA
      }
      NDVI.vals <- c(NDVI.vals,NDVI.val)
    }
  }
  
  fileName <- paste("GOES_NDVI_Diurnal",siteID,"_Oct3_kappaDQF.csv",sep="")
  output <- rbind(t(day.time.vals),NDVI.vals)
  write.table(output,file=fileName,sep=",",col.names=FALSE,row.names=FALSE)
}


siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
iseq <- c(seq(3,6),seq(9,11),seq(15,20))
startDay <- 295
endDay <- 304
output <- 
foreach(i = iseq) %dopar% {
  siteName <- as.character(siteData[i,1])
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
  print(siteName)
  print(c(lat,long))
  createNDVI_GOES_diurnal(lat=lat, long=long, siteID=siteName,startDay=startDay,endDay=endDay)
}


#createNDVI_GOES_diurnal(lat=42.5378, long=-72.1715, siteID="HarvardForest",startDay=366,endDay=250)


