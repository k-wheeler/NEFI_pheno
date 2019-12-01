##' Creates data object with the desired Sf means and variances
##'
##' @param lat The desired latitude if not willowCreek
##' @param long The desired longitude if not willowCreek
##' @param dates The desired days
##' @param siteName The site name
##' @param dataDirectory The data directory
##' @param endDate The end date
##' @param GEFS_Files A vector of the file names for the GEFS files
##' @param GEFS_Directory The directory where the GEFS files are located
##' @param calDatesT
##' @export
createTairs <- function(lat="",long="",dates,siteName,dataDirectory,endDate,GEFS_Files,GEFS_Directory,forecastLength, station,calDatesT=TRUE) {
  years <- lubridate::year(dates)
  if(calDatesT){
    calDates <- seq(dates[1],as.Date(paste((years[length(years)]-1),"-12-31",sep="")),"day")
    ##The sources of the calibration and current measurements differ between willowCreek and other sites
    # if(siteName=="willowcreek"){
    #   TairsCalInd <- download_US_WCr_met(start_date=calDates[1],end_date=calDates[length(calDates)])
    #
    #   ##Need to pad TairsCal to have ensembles
    #   TairsCal <- matrix(ncol=10,nrow=length(TairsCalInd))
    #   for(e in 1:10){
    #     TairsCal[,e] <- TairsCalInd
    #   }
    #   TairsCurrentInd <- download_US_WCr_met(start_date=as.Date("2019-01-01"),end_date=(endDate-forecastLength))
    #
    # }
    #else{
    TairsCal <- load_ERA5_Tair(lat=lat,long=long,years=seq(years[1],2018)) ##columns are each an ensemble (not divided by year)
    print("dim(TairsCal)")
    print(dim(TairsCal))
  }

  TairsCurrentInd <- load_NOAA_met(station=station,startDate=as.Date("2019-08-01"),endDate=(endDate-forecastLength)) ##Array of numeric values

  ##GEFS Forecast (same for all sites) and pad TairsCurrent to be ensembles
  TairsCurrent <- matrix(ncol=length(GEFS_Files),nrow=length(TairsCurrentInd))
  TairsForecast <- numeric()
  for(e in 1:length(GEFS_Files)){
    TairsForecastInd <- load_GEFS_Forecast(dataDirectory=GEFS_Directory,fileName=GEFS_Files[e])

    TairsForecast <- cbind(TairsForecast,TairsForecastInd)
    TairsCurrent[,e] <- TairsCurrentInd
  }
  TairsCurrent <- t(TairsCurrent)
  NOAAmetDays <- seq(as.Date("2019-08-01"),(endDate-forecastLength),"day") ###Changed for autumn
  #print(NOAAmetDays)
  #print(TairsCurrent[,1])
  #print(TairsCurrent[,1]==-9999)
  #print(sum(TairsCurrent[,1]==-9999))
  #print(NOAAmetDays[TairsCurrent[,1]==-9999])
  if(sum(TairsCurrent[,1]==-100)>0){
    #print("filling")
    TairsCurrent[TairsCurrent[,1]< -100,] <- fillNOAAlag(days=NOAAmetDays[TairsCurrent[,1]< -100],siteName=siteName)
    #print("finished filling")
  }
  # if(sum(TairsCurrent[,1]==100)>0){
  #   #print("filling")
  #   TairsCurrent[TairsCurrent[,1]< 100,] <- fillNOAAlag(days=NOAAmetDays[TairsCurrent[,1]< 100],siteName=siteName)
  #   #print("finished filling")
  # }
  ##Create Sfs
  if(calDatesT){
    TairsMeansCal <- rowMeans(TairsCal)
    print('length(TairsMeansCal)')
    print(length(TairsMeansCal))
    TairsVarCal <- apply(TairsCal,MARGIN=1,FUN=var)
  }

  ##Current year
  curDates=seq(as.Date("2019-08-01"),endDate,"day") ##Includes forecasted period ##Changed for autumn

  TairsMeansCur <- c(colMeans(TairsCurrent),rowMeans(TairsForecast))

  TairsVarCur <- c(apply(TairsCurrent,MARGIN=2,FUN=var),apply(TairsForecast,MARGIN=1,FUN=var))

  if(calDatesT){
    TairsMeans <- c(TairsMeansCal,rep(NA,212),TairsMeansCur) ##Pad for the first half of 2019
    TairsVar <- c(TairsVarCal,rep(NA,212),TairsVarCur)
  }
  else{

    TairsMeans <- TairsMeansCur
    TairsVar <- TairsVarCur
  }
  #print(TairsMeans[TairsMeans< (-100)])
  #print(TairsMeans< (-100))
  print(dim(TairsMeans))
  for(t in 1:length(TairsMeans[,ncol(TairsMeans)])){
    if(TairsMeans[t,ncol(TairsMeans)]< (-100)){
      TairsMeans[t,ncol(TairsMeans)] <- mean(TairsMeans[(t-1),ncol(TairsMeans)],TairsMeans[(t+1),ncol(TairsMeans)])
    }
  }

  #TairsVar[is.infinite(TairsVar)] <- 100
  TairsPrec <- 1/TairsVar
  TairsPrec[is.infinite(TairsPrec)] <- 1/100

  #dat <- list(Sf=SfsMeans,Sfprec=1/SfsVar)
  dat <- list(TairMu=TairsMeans,TairPrec=TairsPrec)
  #print("finished createTairs")
  return(dat)
}
