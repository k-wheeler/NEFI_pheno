#' Function to download all MODIS data for input into forecast model
#'
#' @param startDate The start date of data in date format
#' @param endDate The end date of data in date format
#' @param dataDirectory The file path for the directory to download and store data in
#' @param lat Latitide of the site in decimal degrees
#' @param long Longitude of the site in decimal degrees
#' @param siteName The site name for file naming
#' @export
#' @import MODISTools
downloadMODIS <- function(startDate,endDate,dataDirectory,lat,long,siteName){
  DQFfileName <- paste(siteName,"_MOD13Q1_250m_16_days_pixel_reliability_",startDate,endDate,".csv",sep="")
  if(!file.exists(DQFfileName)){
    print("Downloading DQF")
    mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=startDate,end=endDate,site_name = siteName,out_dir = dataDirectory,internal=FALSE)
  }
  NDVIfileName <- paste(siteName,"_MOD13Q1_250m_16_days_NDVI_",startDate,endDate,".csv",sep="")
    #harvard_rel_MOD13Q1_250m_16_days_pixel_reliability_2009-01-012019-09-03.csv
  if(!file.exists(NDVIfileName)){
    print("Downloading NDVI")
    mt_subset(product = "MOD13Q1",lat=lat,lon=long,band=paste("250m_16_days_","NDVI",sep=""),start=startDate,end=endDate,site_name = siteName,out_dir = dataDirectory,internal=FALSE)
  }
  EVIfileName <- paste(siteName,"_MOD13Q1_250m_16_days_EVI_",startDate,endDate,".csv",sep="")
    #harvard_rel_MOD13Q1_250m_16_days_pixel_reliability_2009-01-012019-09-03.csv
  if(!file.exists(EVIfileName)){
    print("Downloading EVI")
    mt_subset(product = "MOD13Q1",lat=lat,lon=long,band=paste("250m_16_days_","EVI",sep=""),start=startDate,end=endDate,site_name = siteName,out_dir = dataDirectory,internal=FALSE)
  }

  #newDQFFileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new DQF data downloaded
  #print(newDQFFileName)
  newDQFDat <- read.csv(DQFfileName,header=TRUE) ##The new DQF data
  DQFdata <- newDQFDat$value

  newDat <- read.csv(NDVIfileName,header=TRUE)
  newNDVIdat <- cbind(newDat,DQFdata)

  newDat <- read.csv(EVIfileName,header=TRUE)
  newEVIdat <- cbind(newDat,DQFdata)

  fileName <- paste(dataDirectory,siteName,"_",NDVI,"_MOD13Q1_",startDate,"_",endDate,".csv",sep="")
  write.table(newNDVIdat,file=fileName,row.names = FALSE,col.names = TRUE,sep=",")

  fileName <- paste(dataDirectory,siteName,"_",EVI,"_MOD13Q1_",startDate,"_",endDate,".csv",sep="")
  write.table(newEVIdat,file=fileName,row.names = FALSE,col.names = TRUE,sep=",")
}
# write.table(dat,file=fileName,row.names = FALSE,col.names = TRUE,sep=",")
#
#
#
#
#
#
#
#   fileName <- paste(dataDirectory,siteName,"_",metric,"_MOD13Q1_",startDate,"_",endDate,".csv",sep="")
#   print(fileName)
#   if(!file.exists(fileName)){
#     #files <- intersect(intersect(dir(path=dataDirectory,pattern=paste(siteName,"_",metric,sep="")),dir(path=dataDirectory,pattern="MOD13Q1")),dir(path=dataDirectory,pattern=as.character(startDate))) #Current downloaded data files
#     #if(length(files>0)){ #If there is a data file, identify what the last date you have downloaded
#     #  firstDate <- strsplit(files[length(files)],split="_")[[1]][4]
#     #  lastDate <- strsplit(strsplit(files[length(files)],split="_")[[1]][5],split=".c")[[1]][1]
#     #}
#     #else{
#     lastDate <- (as.Date(startDate) - 1) #If no data files have been downloaded, the last date you have downloaded is your start date - 1
#     #}
#     print(lastDate)
#     # if(metric!="rel"){
#     #   newDQFFileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new DQF data downloaded
#     #   if(!file.exists(newDQFFileName)){
#     #     print("Downloading MODIS DQF File")
#     #     try(mt_subset(product = "MOD13Q1",lat=lat,lon=long,
#     #                   band="250m_16_days_pixel_reliability",start=(lastDate+1),end=endDate,
#     #                   site_name = paste(siteName,"_rel",sep=""),out_dir = dataDirectory,internal=FALSE),silent=TRUE)
#     #   }
#     # }
#     print("Downloading MODIS Index File")
#     if(metric =="rel"){
#       #mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=(lastDate+1),end=endDate,site_name = paste(siteName,"_",metric,sep=""),out_dir = dataDirectory,internal=FALSE)
#     }
#     else{
#       #mt_subset(product = "MOD13Q1",lat=lat,lon=long,band=paste("250m_16_days_",metric,sep=""),start=(lastDate+1),end=endDate,site_name = paste(siteName,"_",metric,sep=""),out_dir = dataDirectory,internal=FALSE)
#     }
#     newFileName <- paste(dataDirectory,siteName,"_",metric,"_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new data downloaded
#
#     #if(length(files>0)){
#     #  dat <- read.csv(paste(dataDirectory,files[length(files)],sep=""),header=TRUE) ##Reads the old data file
#     #  if(file.exists(newFileName)){ ##If additional available MODIS data was downloaded
#     #    newDat <- read.csv(newFileName,header=TRUE,skip=15) ##The new data
#     #    if(metric!="rel"){
#     #      newDQFDat <- read.csv(newDQFFileName,header=TRUE,skip=15) ##The new DQF data
#     #      DQFdata <- newDQFDat$data
#     #      newDat <- cbind(newDat,DQFdata)
#     #    }
#     #    dat <- rbind(dat,newDat) ##Combine the new and old data files
#     #  }
#     #}
#     # else{ #If there were no data files to start with than your new data file will be your final one
#     newDat <- read.csv(newFileName,header=TRUE,skip=15)
#     #print(newDat)
#     if(metric!="rel"){
#       newDQFFileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",(as.Date(lastDate)+1),"_",endDate,".csv",sep="") #File name for new DQF data downloaded
#       print(newDQFFileName)
#       newDQFDat <- read.csv(newDQFFileName,header=TRUE) ##The new DQF data
#       print(newDQFDat$data)
#       DQFdata <- newDQFDat$data
#       newDat <- cbind(newDat,DQFdata)
#     }
#     dat <- newDat
#     #}
#     write.table(dat,file=fileName,row.names = FALSE,col.names = TRUE,sep=",")
#   }
# }
