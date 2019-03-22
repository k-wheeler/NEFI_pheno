#' Function to download MODIS data for input into forecast model
#'
#' @param startDate The start date of data in date format
#' @param endDate The end date of data in date format
#' @param metric Desired metric to download (NDVI or EVI)
#' @param dataDirectory The file path for the directory to download and store data in
#' @param lat Latitide of the site in decimal degrees
#' @param long Longitude of the site in decimal degrees
#' @param siteName The site name for file naming
#' @export
#' @import MODISTools
downloadMODIS <- function(startDate,endDate,metric,dataDirectory,lat,long,siteName){
  fileName <- paste(dataDirectory,siteName,"_",metric,"_MOD13Q1_",startDate,"_",endDate,".csv",sep="")
  if(!file.exists(fileName)){
    files <- intersect(intersect(dir(path=dataDirectory,pattern=paste(siteName,"_",metric,sep="")),dir(path=dataDirectory,pattern="MOD13Q1")),dir(path=dataDirectory,pattern=startDate)) #Current downloaded data files
    if(length(files>0)){ #If there is a data file, identify what the last date you have downloaded
      firstDate <- strsplit(files[length(files)],split="_")[[1]][4]
      lastDate <- strsplit(strsplit(files[length(files)],split="_")[[1]][5],split=".c")[[1]][1]
    }
    else{
      lastDate <- (as.Date(startDate) - 1) #If no data files have been downloaded, the last date you have downloaded is your start date - 1
    }
    newDQFFileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new DQF data downloaded
    if(!file.exists(newDQFFileName)){
      print("Downloading MODIS DQF File")
      try(mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=(lastDate+1),end=endDate,site_name = paste(siteName,"_rel",sep=""),out_dir = dataDirectory,internal=FALSE),silent=TRUE)
    }
     print("Downloading MODIS Index File")
    try(mt_subset(product = "MOD13Q1",lat=lat,lon=long,band=paste("250m_16_days_",metric,sep=""),start=(lastDate+1),end=endDate,site_name = paste(siteName,"_",metric,sep=""),out_dir = dataDirectory,internal=FALSE),silent=TRUE)
    newFileName <- paste(dataDirectory,siteName,"_",metric,"_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new data downloaded

    if(length(files>0)){
      dat <- read.csv(paste(dataDirectory,files[length(files)],sep=""),header=TRUE) ##Reads the old data file
      if(file.exists(newFileName)){ ##If additional available MODIS data was downloaded
        newDat <- read.csv(newFileName,header=TRUE,skip=15) ##The new data
        newDQFDat <- read.csv(newDQFFileName,header=TRUE,skip=15) ##The new DQF data
        DQFdata <- newDQFDat$data
        newDat <- cbind(newDat,DQFdata)
        dat <- rbind(dat,newDat) ##Combine the new and old data files
      }
    }
    else{ #If there were no data files to start with than your new data file will be your final one
      newDat <- read.csv(newFileName,header=TRUE,skip=15)
      newDQFDat <- read.csv(newDQFFileName,header=TRUE,skip=15) ##The new DQF data
      DQFdata <- newDQFDat$data
      newDat <- cbind(newDat,DQFdata)
      dat <- newDat
    }
    write.table(dat,file=fileName,row.names = FALSE,col.names = TRUE,sep=",")
  }
}
