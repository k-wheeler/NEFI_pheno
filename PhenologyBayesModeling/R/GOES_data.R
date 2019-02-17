##' For GOES NDVI data, construct the data object for input into MCMC
##'
##' @param siteName Site Name
##' @param lat Latitude
##' @param long Longitude
##' @param startDay The start day counted as the day number after 2016-12-31
##' @param endDay The end day counted as the day number after 2016-12-31
##' @param TZ The timezone off of UTC
##' @param window Boolean of if you want to average out over a 4 hour window
##' @export
GOES_data <- function(siteName,lat,long,startDay,endDay,TZ,window=FALSE,maxValue=FALSE) {
  startDate <- as.Date(startDay,origin="2016-12-31")
  endDate <- as.Date(endDay,origin="2016-12-31")
  if(window){
    fileName <- paste("GOES_NDVI_",siteName,"_",startDate,"_",endDate,"_Avg.csv",sep="")
    #fileName <- "GOES_NDVI_HarvardForest_Avg.csv"
    if(!file.exists(fileName)){
      createNDVI_GOES_Avg(lat=lat,long=long,startDay=startDay,endDay=endDay,fileName=fileName,TZ=TZ)
    }
  }
  else if(maxValue){
    fileName <- paste("GOES_NDVI_",siteName,"_",startDate,"_",endDate,"_max.csv",sep="")
  }
  else{
    fileName <- paste("GOES_NDVI_",siteName,"_",startDate,"_",endDate,"_noon.csv",sep="")
    if(!file.exists(fileName)){
      createNDVI_GOES(lat=lat,long=long,startDay=startDay,endDay=endDay,fileName=fileName,TZ=TZ)
    }
  }
  print(fileName)
  GOES <- read.csv(fileName,header=FALSE)
  print(dim(GOES))
  GOES <- GOES[,colSums(is.na(GOES)) == 0]
  GOES_Days <- as.numeric(GOES[1,])
  GOES_NDVI <- as.numeric(GOES[2,])
  print(GOES_Days[1:20])

  # for(i in 1:length(GOES_Days)){
  #   if(GOES_Days[i]<100){
  #     GOES_Days[i] <- as.numeric(GOES_Days[i]) + 365
  #   }
  # }
  # bk <- which(GOES_Days>=366)[1]
  # print(bk)
  # for(i in bk:length(GOES_Days)){
  #   GOES_Days[i] <- GOES_Days[i]+365
  # }

  #days <- which(GOES_Days>startDay&GOES_Days<endDay)

  y <- GOES_NDVI
  x <- GOES_Days

  data <- list(x=x,y=y,n=length(y))
  if(window){
    data$obs.prec <- as.numeric(GOES[3,])
    print(length(as.numeric(GOES[3,])))
    #print(data$obs.prec)
    data$size <- as.numeric(GOES[4,])
  }
  return(data)
}

