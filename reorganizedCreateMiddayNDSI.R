#!/usr/bin/env Rscript

install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
library("ncdf4")
library(plyr)
library("PhenologyBayesModeling")
library(doParallel)

#detect cores.
n.cores <- 6

#register the cores.
registerDoParallel(cores=n.cores)

##Can only do for same TZ

createNDSI_sub <- function(siteData,orbitVersion,day.time,year){
  print(orbitVersion)
  ##Will need to return a vector of the NDSI values
  
  ##Create File Paths
  ACM.path <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=paste("OR_ABI-L2-ACMC-M3_G16_s",day.time,sep="")),sep="")
  filestrC05 <- paste("OR_ABI-L1b-RadC-M3C05_G16_s",day.time,sep="")
  #print(filestrC05)
  filestrC02 <- paste("OR_ABI-L1b-RadC-M3C02_G16_s",day.time,sep="")
  filePathC02 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC02),sep="")
  filePathC05 <- paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrC05),sep="")
  NDSIs <- numeric()
  #print(ACM.path)
  #print(filePathC02)
  #print(filePathC05)
  
  ##Open Files
  #print(paste("filePathC02:",filePathC02))
  #print(paste("filePathC05:",filePathC05))
  if(nchar(filePathC02)>20 && nchar(filePathC05)>20){
    if(file.exists(filePathC02) && file.exists(filePathC05) && !dir.exists(filePathC02) && !dir.exists(filePathC05)){
      ACM.file <-nc_open(ACM.path)
      R2.file <- nc_open(filePathC02)
      R5.file <- nc_open(filePathC05)
      
      ##Extract Data
      R5 <- ncvar_get(R5.file,"Rad")
      R2 <- ncvar_get(R2.file,"Rad") #the full R2 dataset
      R5.kappa0 <- ncvar_get(R5.file,"kappa0")
      R2.kappa0 <- ncvar_get(R2.file,"kappa0")
      R5.DQF <- ncvar_get(R5.file,"DQF") #Data Quality Flags
      R2.DQF <- ncvar_get(R2.file,"DQF")
      R5 <- R5 * R5.kappa0 #done to covert radiance to reflectance
      R2 <- R2 * R2.kappa0
      clouds <- ncvar_get(ACM.file,"BCM")
      clouds.DQF <- ncvar_get(ACM.file,"DQF")
      print(dim(R5.DQF))
      
      for(i in 1:nrow(siteData)){
        ##General Site Data
        siteName <- as.character(siteData[i,1])
        lat <- as.numeric(siteData[i,2])
        long <- as.numeric(siteData[i,5])
        TZ <- as.numeric(siteData[i,6])
        
        ##Determine index values
        lat.rd <- lat*2*pi/360
        long.rd <- long*2*pi/360
        
        Ind2 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion=orbitVersion),2,orbitVersion=orbitVersion)
        Ind5 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion=orbitVersion),5,orbitVersion=orbitVersion)
        ACM.ind <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion=orbitVersion),"ACM",orbitVersion=orbitVersion)
        print(paste("Ind5:",Ind5))
        i2 <- Ind2[1]
        j2 <- Ind2[2]
        i5 <- Ind5[1]
        j5 <- Ind5[2]
        
        if(!is.na(R5.DQF[i5,j5]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[i2,j2]) && !is.na(R2.DQF[(i2+1),j2]) && !is.na(R2.DQF[i2,(j2+1)]) && !is.na(R2.DQF[(i2+1),(j2+1)]) && !is.na(clouds[ACM.ind[1],ACM.ind[2]]) && !is.na(clouds.DQF[ACM.ind[1],ACM.ind[2]])){
          if(R5.DQF[i5,j5]==0 && R2.DQF[i2,j2]==0 && R2.DQF[i2,j2]==0 && R2.DQF[(i2+1),j2]==0 && R2.DQF[i2,(j2+1)]==0 && R2.DQF[(i2+1),(j2+1)]==0 && clouds[ACM.ind[1],ACM.ind[2]] == 0 && clouds.DQF[ACM.ind[1],ACM.ind[2]] == 0){
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
        print(paste("NDSI:",output))
        NDSIs <- c(NDSIs,output)
      }
    }
    else{
      return(rep(NA,nrow(siteData)))
    }
  }
  else{
    return(rep(NA,nrow(siteData)))
  }
  return(NDSIs)
}

createEmptyFiles <- function(siteData,day,year){
  for(i in 1:length(siteData)){
    siteName <- as.character(siteData[i,1])
    fileName <- paste("GOES_NDSI_Midday",siteName,"_",day,"_",day,"_kappaDQF.csv",sep="")
    write.table(NA,file=fileName,sep=",",col.names=FALSE,row.names=FALSE)
  }
}

checkFileExists <- function(siteName,day){
  ##Returns True if the file exists
  dat <- read.csv(file=paste(siteName,"_missingMidday_NDSI_Days.csv",sep=""),header=FALSE)
  for(i in 1:nrow(dat)){
    if(as.numeric(dat[i,1])==as.numeric(day)){
      return(FALSE)
    }
  }
  return(TRUE)
}

createMissingFilesList <- function(siteName){
  diurnal.files <- dir(pattern=paste("GOES_NDSIMidday",siteName,sep=""))
  days <- numeric()
  #print(paste("Length of diurnal.files:"),length(diurnal.files))
  if(length(diurnal.files)!= 0){
    for(i in 1:length(diurnal.files)){
      #print(diurnal.files[i])
      start <- as.numeric(strsplit(diurnal.files[i],"_")[[1]][4])
      end <- as.numeric(strsplit(diurnal.files[i],"_")[[1]][5])
      #print(start)
      #print(end)
      dys <- seq(start,end,1)
      days <- c(days,dys)
    }
  }
  #sort(days)
  #all.days <- c(seq(321,333,1),seq(348,365,1))
  all.days <- c(seq(1,333,1),seq(348,365,1))
  #all.days <- seq(1,181,1)
  missingDays <- numeric()
  for(j in 1:length(all.days)){
    if((!all.days[j] %in% days)){
      missingDays <- c(missingDays,all.days[j])
    }
  }
  #missingDays
  #print(length(missingDays))
  if(length(missingDays)==0){
    missingDays <- c(missingDays,0)
  }
  write.table(missingDays,file=paste(siteName,"_missingMidday_NDSI_Days.csv",sep=""),sep=",",col.names=FALSE,row.names=FALSE)
}


createNDSI_GOES_MiddayMAIN <- function(day,siteData,orbitVersion,year,TZ){
  print(orbitVersion)
  #print("Inside Main")
  ##hours for Midday calculations
  #print(TZ)
  hrs <- as.character(seq(10+as.numeric(TZ),(13+as.numeric(TZ)),1))
  #print(hrs)
  if(as.numeric(day)<10){
    day <- paste("00",as.character(day),sep="")
  }
  else if(as.numeric(day)<100){
    day <- paste("0",as.character(day),sep="")
  }
  else{
    day <- as.character(day)
  }
  ##
  filestrACM <- paste("OR_ABI-L2-ACMC-M3_G16_s",year,day,sep="")
  #print(filestrACM)
  ACM.files <- character()
  for(q in 1:length(hrs)){
    if(as.numeric(hrs[q])<10){
      hrs[q] <- paste("0",hrs[q],sep="")
    }
    else if(as.numeric(hrs[q])>23){
      hrs[q] <- paste("0",as.character((hrs[q]-24)),sep="")
      day <- as.numeric(day)+1
      if(as.numeric(day)<10){
        day <- paste("00",as.character(day),sep="")
      }
      else if(as.numeric(day)<100){
        day <- paste("0",as.character(day),sep="")
      }
      else{
        day <- as.character(day)
      }
    }
    #print(c("hrs[q]",hrs[q]))
    newFiles <- intersect(dir(path="GOES_Data2017",pattern=filestrACM),dir(path="GOES_Data2017",pattern=paste("s",year,day,hrs[q],sep="")))
    ACM.files <- c(ACM.files,newFiles)
  }
  #print(ACM.files)
  
  if(!dir.exists((paste("GOES_Data2017/",dir(path="GOES_Data2017",pattern=filestrACM),sep="")))){
    if(length(ACM.files>1)){
      day.time.vals <- character()
      NDSI.vals <- matrix(ncol=nrow(siteData),nrow=0)
      #print(dim(NDSI.vals))
      for(j in 1:length(ACM.files)){
        ##Open Files
        day.time <- substr(ACM.files[j],24,34)
        day.time.vals <- c(day.time.vals,day.time)
        #print(NDSI)
        NDSI.vals <- rbind(NDSI.vals,createNDSI_sub(siteData=siteData,orbitVersion=orbitVersion,day.time=day.time) )
        #print(NDSI.vals)
      }
      print(dim(NDSI.vals))
      for(i in 1:nrow(siteData)){
        siteName <- as.character(siteData[i,1])
        fileName <- paste("GOES_NDSI_Midday",siteName,"_",day,"_",day,"_kappaDQF.csv",sep="")
        #print("Printing NDSI.vals[,i]:")
        #print(NDSI.vals[,i])
        #print(dim(NDSI.vals[,i]))
        output <- rbind(day.time.vals,NDSI.vals[,i])
        write.table(output,file=fileName,sep=",",col.names=FALSE,row.names=FALSE)
      }
    }
    else{
      createEmptyFiles(siteData=siteData,day=day)
    }
  }
  else{
    createEmptyFiles(siteData=siteData,day=day)
  }
}

#########

##For TZ == 5 
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)[c(1,2,3,4,5,8,11,16,17,18,19,20,21,22),] ##TZ 5 sites
#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)[c(8,11,16,17,18,19,20),]#[c(2,3,4),]
#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)[c(1,2,3,4,5),]
#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)[c(6,9,10,15),] ##TZ 6 sites
#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)[c(7,13,14),]##TZ 7 sites
#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)[c(12),] ##TZ 8
#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)[c(21,22),]
TZ <- as.numeric(siteData[1,6])
PFT <- as.character(siteData[1,5])
for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData[s,1])
  #print(siteName)
  createMissingFilesList(siteName)
}
print("Done creating missing day files")
#all.days <- c(seq(321,333,1),seq(348,365,1))
all.days <- c(seq(1,151,1),seq(274,333,1),seq(348,365,1))
#all.days <- seq(158,181)
#year <- 2017
#output <- 
#foreach (d = 1:length(all.days)) %dopar% {
for(d in 1:length(all.days)){
  print(paste("Starting Day:",all.days[d],sep=" "))
  iseq <- numeric()
  if(PFT=="SH"){
    if(all.days[d]<110){#######Need to change for DB
      year <- 2018
    }
    else{
      year <- 2017
    }
    if(all.days[d]<334 && all.days[d]>109){##Need to change for DB
      orbitVersion <- "OLD"
    }
    else{
      orbitVersion <- "NEW"
    }
  }
  else if(PFT=="DB"){
    if(all.days[d]<182){
      year <- 2018
    }
    else{
      year <- 2017
    }
    if(all.days[d]<334 && all.days[d]>181){
      orbitVersion <- "OLD"
    }
    else{
      orbitVersion <- "NEW"
    }
  }
  for(s in 1:nrow(siteData)){
    siteName <- as.character(siteData[s,1])
    #print(siteName)
    if(!checkFileExists(siteName = siteName,day=all.days[d])){
      iseq <- c(iseq,s)
    }
  }
  print(iseq)

  if(length(iseq)>0){
    createNDSI_GOES_MiddayMAIN(day=all.days[d],siteData=siteData[iseq,],orbitVersion = orbitVersion,year = year,TZ = TZ)
  }
  
}

