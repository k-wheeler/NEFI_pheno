##' For MODIS EVI data, construct the data object for input into MCMC
##'
##' @param siteName Site Name
##' @param lat Latitude
##' @param long Longitude
##' @param startDate The start day counted as the day number after 2016-12-31
##' @param endDate The end day counted as the day number after 2016-12-31
##' @param metric "NDVI" or "EVI"
##' @param dataDirectory The data directory
##' @import MODISTools
##' @export
MODIS_data <- function(siteName,lat,long,startDate="",endDate="",metric,startDay=FALSE,endDay=FALSE,lastYear=2018,dataDirectory="") {
  if(typeof(startDay)!=typeof(FALSE)){
    startDate <- as.Date(startDay,origin="2016-12-31")
    endDate <- as.Date(endDay,origin="2016-12-31")
  }

  fileName <- paste(dataDirectory,siteName,"_",metric,"_MOD13Q1_",startDate,"_",endDate,".csv",sep="")
  print(fileName)
  # if(!file.exists(fileName)){
  #   print("Downloading MODIS File")
  #   directory=getwd()
  #   mt_subset(product = "MOD13Q1",lat=lat,lon=long,band=paste("250m_16_days_",metric,sep=""),start=startDate,end=endDate,site_name = paste(siteName,"_",metric,sep=""),out_dir = directory,internal=FALSE)
  # }
  #DQFfileName <- paste(dataDirectory,siteName,"_","rel","_MOD13Q1_",startDate,"_",endDate,".csv",sep="")
  # if(!file.exists(DQFfileName)){
  #   print("Downloading MODIS DQF File")
  #   directory=getwd()
  #   mt_subset(product = "MOD13Q1",lat=lat,lon=long,band="250m_16_days_pixel_reliability",start=startDate,end=endDate,site_name = paste(siteName,"_","rel",sep=""),out_dir = directory,internal=FALSE)
  # }
  # print("MODIS File Downloaded")
  # dat <- read.csv(fileName,header=TRUE,skip=15)
  # print(dat)
  #DQF <- read.csv(DQFfileName,header=TRUE)
  dat <- read.csv(fileName,header=TRUE)
  x <- numeric()
  y <- numeric()

  for(i in 1:nrow(dat)){
    #print(i)
    y <- c(y,dat$value[i]/10000)
    #DQF.val <- DQF$data[i]
    DQF.val <- dat$DQFdata[i]
    #print(DQF.val)
    if(DQF.val!= 0 && DQF.val != 1){
      y[i] <- NA
    }
    #print(DQF.val)
    #print(y[i])
    #print(dat$calendar_date[i])
    tmp <- as.Date(dat$calendar_date[i])
    #print(tmp)
    # #yr <- strsplit(preTmp,"/")[[1]][3]
    # yr <- lubridate::year(as.Date(preTmp))
    # print(yr)
    # mth <- strsplit(preTmp,"/")[[1]][1]
    # mth <- lubridate::month(as.Date(preTmp))
    # #print(mth)
    # #dy <- strsplit(preTmp,"/")[[1]][2]
    # dy <- lubridate::day(as.Date(preTmp))
    # print(dy)
    # if(as.numeric(mth)<10){
    #   mth <- paste("0",mth,sep="")
    # }
    # if(as.numeric(dy)<10){
    #   dy <- paste("0",dy,sep="")
    # }
    # preTmpDate <- paste("20",yr,"-",mth,"-",dy,sep="")
    # tmp <- as.Date(preTmpDate)
    #print(tmp)
    x.val <- as.numeric(format(tmp, "%j"))
    if(substr(tmp,1,4)==lastYear){
      x.val <- x.val + 365
    }
    x <- c(x,x.val)
  }

  data <- list(x=x,y=y,n=length(y))
  return(data)
}
