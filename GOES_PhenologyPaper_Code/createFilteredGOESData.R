#' checkHours
#'
#' @param dayData 
#'
#' @return
#' @export
#'
#' @examples
checkHours <- function(dayData){
  hours <- seq(4,24)
  checks <- 0
  for(h in hours){
    hourData <- as.matrix(dayData[,floor(dayData[3,])==h])
    if(ncol(hourData)>0){
      if(length(na.omit(as.numeric(hourData[2,])))>0){
        checks <- checks + 1
      }
    }
  }
  return(checks)
}
#' filterData
#'
#' @param dayData 
#' @param date.val 
#' @param lat 
#' @param long 
#' @param TZ_name 
#'
#' @return
#' @export
#'
#' @examples
filterData <- function(dayData,date.val,lat,long,TZ_name){
  suntimes <- getSunlightTimes(date=date.val,lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),tz = TZ_name)
  dawnTime <- lubridate::hour(suntimes$nauticalDawn)+(lubridate::minute(suntimes$nauticalDawn)/60)
  duskTime <- lubridate::hour(suntimes$nauticalDusk)+(lubridate::minute(suntimes$nauticalDusk)/60)
  
  for(c in 1:ncol(dayData)){
    if(!is.na(dayData[2,c])){
      if(round(dayData[2,c],digits=4)==0.6040){
        dayData[2,c] <- NA
      }else if(dayData[2,c]>0.9){
        dayData[2,c] <- NA
      }else if(dayData[2,c]<0){
        dayData[2,c] <- NA
      }
    }
  }
  
  dayData[2,as.numeric(dayData[3,])<(as.numeric(dawnTime)+1.5)] <- NA
  dayData[2,as.numeric(dayData[3,])>(as.numeric(duskTime)-1.5)] <- NA
  return(dayData)
}


#' createFilteredData
#'
#' @param DiurnalFitSavePath
#' @param siteData
#' @param savePath
#' @param calculatedNDVIGOESpath
#' @param yr
#' @param startDate
#' @param endDate
#' 
#' @return
#' @import ncdf4
#' @import plyr
#' @import PhenologyBayesModeling
#' @import rjags
#' @import runjags
#' @export
#'
#' @examples
createFilteredGOESData <- function(DiurnalFitSavePath,siteData,savePath,calculatedNDVIGOESpath,yr,startDate,endDate){
  # library(doParallel)
  # n.cores <- 15
  # registerDoParallel(cores=n.cores)
  # output <- 
  #   foreach(s = iseq) %dopar% {
  for(s in 1:nrow(siteData)){
    siteName <- as.character(siteData[s,1])
    print(siteName)
    lat <- as.numeric(siteData$Lat[i])
    long <- as.numeric(siteData$Long[i])
    TZ <- as.character(siteData$TZ[i])
    if(TZ==5){
      TZ_name <- "America/New_York"
    }else if(TZ==6){
      TZ_name <- "America/Chicago"
    }else if(TZ_name==""){
      print("TZ_name unknown. Please enter it")
      return(FALSE)
    }
    outDataFile <- paste(savePath,siteName,"_",yr,"_diurnalFitDataFiltered.RData",sep="")

    diurnalFits <- intersect(dir(path=DiurnalFitSavePath,pattern="varBurn"),
                             dir(path=DiurnalFitSavePath,
                                 pattern=paste(siteName,"_",yr,sep="")))
    c.vals <- numeric()
    prec.vals <- numeric()
    days <- numeric()
    ##Counts total number of observations for each day/site
    cts <- numeric()
    for(i in 1:length(diurnalFits)){
      print(diurnalFits[i])
      load(paste(DiurnalFitsavePath,diurnalFits[i],sep=""))
      if(typeof(var.burn)!=typeof(FALSE)){
        out.mat <- data.frame(as.matrix(var.burn))
        print(colnames(out.mat))
        c <- mean(out.mat$c)
        prec <- 1/var(out.mat$c)
        
        dy <- strsplit(diurnalFits[i],"_")[[1]][2]
        
        dayDataFile <- intersect(dir(path=calculatedNDVIGOESpath,pattern=siteName),
                                 dir(path=calculatedNDVIGOESpath,pattern=paste(dy,".csv",sep="")))
        print(dayDataFile)
        dayData <- read.csv(paste(calculatedNDVIGOESpath,dayDataFile,sep=""),header=FALSE)
        ct <- length(dayData[2,][!is.na(dayData[2,])])
        if(ct>1){
          c.vals <- c(c.vals,c)
          prec.vals <- c(prec.vals,prec)
          days <- c(days,substr(dy,5,7))
          cts <- c(cts,ct)
        }
      }
    }
    save(cts, file=paste(savePath,siteName,"_",yr,"_counts.RData",sep=""))
    data <- list()
    data$x <- as.numeric(days)
    data$y <- as.numeric(c.vals)
    data$obs.prec <- as.numeric(prec.vals)
    data$n <- length(data$x)
    
    data$Q1 <- data$y - 1.96 * sqrt(1/data$obs.prec)
    data$Q3 <- data$y + 1.96 * sqrt(1/data$obs.prec)
    ctoff <- 24 ###Keep days with more than 24 observations
    enoughHours <- logical()
    
    ##Filters days with not enough hours that have an observation
    for(d in data$x){
      if(as.numeric(d)<10){
        dy <- paste("00",as.character(d),sep="")
      }else if(as.numeric(d)<100){
        dy <- paste("0",d,sep="")
      }else{
        dy <- d
      }
      dayDataFile <- intersect(dir(path=calculatedNDVIGOESpath,pattern=siteName),
                               dir(path=calculatedNDVIGOESpath,pattern=paste(yr,dy,".csv",sep="")))
      dayData <- read.csv(paste(calculatedNDVIGOESpath,dayDataFile,sep=""),header=FALSE)
      date.val <- as.Date(as.numeric(dy),origin=as.Date(paste(as.character(as.numeric(yr)-1),"-12-31",sep="")))
      dayData <- filterData(date.val = date.val,lat=lat,long=long,TZ_name = TZ_name, dayData=dayData)
      hourChecks <- checkHours(dayData = dayData)
      if(hourChecks>4){
        enoughHours <- c(enoughHours,TRUE)
      }else{
        enoughHours <- c(enoughHours,FALSE)
      }
    }
    dataMat <- data.frame(cbind(data$x,data$y,data$obs.prec,data$Q1,data$Q3,cts=cts))
    colnames(dataMat) <- c("x","y","obs.prec","Q1","Q3","cts")
    dataMat <- data.frame(dataMat[enoughHours,])
    dataMat <- data.frame(dataMat[dataMat$cts>ctoff,]) #Removes days that do not have enough observations
    allCounts <- allCounts + length(dataMat$x)
    
    if(yr==2019){ #Removes the problematic days for specific sites
      if(i==1){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(17,306),])
      }else if(i==3){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(51,259,337),])
      }else if(i==4){
        dataMat <- data.frame(dataMat[!dataMat$x==335,])
      }else if(i==11){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(52,62),])
      }else if(i==15){
        dataMat <- data.frame(dataMat[!dataMat$x==11,])
      }else if(i==10){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(52,127,284,286,345,364),])
      }else if(i==14){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(85,82,327),])
      }
    }else if(yr==2018){
      if(i==1){
        dataMat <- data.frame(dataMat[!dataMat$x==78,])
      }else if(i==4){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(71,324,352),])
      }else if(i==7){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(28,359),])
      }else if(i == 9){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(80,90,113,114,123,125),])
      }else if(i==12){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(69,112),])
      }else if(i==13){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(2,74,85,358),])
      }else if(i==15){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(120,124),])
      }else if(i==3){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(43,45,47,61,73,7578,100,101,110,111,112,124,126,139,152,278,284),])
      }else if(i==10){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(38,80,99,108,109,168,208),])
      }else if(i==14){
        dataMat <- data.frame(dataMat[!dataMat$x%in%c(91,97,99,112,113,114,124),])
      }
    }
    
    filteredData <- list(x=dataMat$x,y=dataMat$y,obs.prec=dataMat$obs.prec,
                         Q1=dataMat$Q1,Q3=dataMat$Q3,n=nrow(dataMat))
    
    filteredData <- list(x=dataMat$x,y=dataMat$y,obs.prec=dataMat$obs.prec,n=nrow(dataMat))

    save(filteredData,file=outDataFile)
    
    print("Done with creating Data")
  }
}
