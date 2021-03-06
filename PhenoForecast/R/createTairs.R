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
##' @param TZ_name
##' @param calDatesT
##' @export
createTairs <- function(lat="",long="",dates,siteName,dataDirectory,endDate,GEFS_Files,GEFS_Directory,forecastLength,station,calDatesT=TRUE,TZ_name="") {
  ERA5dataFolder=paste("/projectnb/dietzelab/kiwheel/ERA5/Data/",siteName,"/",sep="")
  years <- lubridate::year(dates)
  ##Calibration
  if(calDatesT){ #Include calibration dates
    TairsCal <- load_ERA5_Tair_New(lat=lat,long=long,endDate=endDate,calDatesT=TRUE,ERA5dataFolder,TZ_name=TZ_name)#Returns a matrix with 10 ensembles each as a row
    print("dim(TairsCal)")
    print(dim(TairsCal))
  }
  TairsMeans <- apply(TairsCal,MARGIN = 2,mean)
  TairsPrec <- 1/apply(TairsCal,MARGIN = 2,var)

  # TairsCurrentInd <- load_NOAA_met(station=station,startDate=as.Date("2019-08-01"),endDate=(endDate-forecastLength)) ##Array of numeric values
  #
  # ##GEFS Forecast (same for all sites) and pad TairsCurrent to be ensembles
  # TairsCurrent <- matrix(ncol=length(GEFS_Files),nrow=length(TairsCurrentInd))
  # TairsForecast <- numeric()
  # for(e in 1:length(GEFS_Files)){
  #   TairsForecastInd <- load_GEFS_Forecast(dataDirectory=GEFS_Directory,fileName=GEFS_Files[e])
  #
  #   TairsForecast <- cbind(TairsForecast,TairsForecastInd)
  #   TairsCurrent[,e] <- TairsCurrentInd
  # }
  # TairsCurrent <- t(TairsCurrent)
  # NOAAmetDays <- seq(as.Date("2019-08-01"),(endDate-forecastLength),"day") ###Changed for autumn
  # #print(NOAAmetDays)
  # #print(TairsCurrent[,1])
  # #print(TairsCurrent[,1]==-9999)
  # #print(sum(TairsCurrent[,1]==-9999))
  # #print(NOAAmetDays[TairsCurrent[,1]==-9999])
  # if(sum(TairsCurrent[,1]==-100)>0){
  #
  #   TairsCurrent[TairsCurrent[,1]< -100,] <- fillNOAAlag(days=NOAAmetDays[TairsCurrent[,1]< -100],siteName=siteName)
  #
  # }
  # # if(sum(TairsCurrent[,1]==100)>0){
  # #   #print("filling")
  # #   TairsCurrent[TairsCurrent[,1]< 100,] <- fillNOAAlag(days=NOAAmetDays[TairsCurrent[,1]< 100],siteName=siteName)
  # #   #print("finished filling")
  # # }
  # ##Create Sfs
  # if(calDatesT){
  #   TairsMeansCal <- rowMeans(TairsCal)
  #   print('length(TairsMeansCal)')
  #   print(length(TairsMeansCal))
  #   TairsVarCal <- apply(TairsCal,MARGIN=1,FUN=var)
  # }
  #
  # ##Current year
  # curDates=seq(as.Date("2019-08-01"),endDate,"day") ##Includes forecasted period ##Changed for autumn
  #
  # TairsMeansCur <- c(colMeans(TairsCurrent),rowMeans(TairsForecast))
  #
  # TairsVarCur <- c(apply(TairsCurrent,MARGIN=2,FUN=var),apply(TairsForecast,MARGIN=1,FUN=var))
  # if(calDatesT){
  #   print("length(TairsMeansCur)")
  #   print(length(TairsMeansCur))
  #   TairsMeans <- c(TairsMeansCal,rep(NA,212),TairsMeansCur) ##Pad for the first half of 2019
  #   TairsVar <- c(TairsVarCal,rep(NA,212),TairsVarCur)
  # }
  # else{
  #   TairsMeans <- TairsMeansCur
  #   TairsVar <- TairsVarCur
  # }
  # #print(TairsMeans[TairsMeans< (-100)])
  # #print(TairsMeans< (-100))
  #
  # for(t in 1:length(TairsMeans)){
  #   #print(t)
  #   #print(TairsMeans[t])
  #   if(!is.na(TairsMeans[t])){
  #     if(TairsMeans[t]< (-100)){
  #       TairsMeans[t] <- mean(TairsMeans[(t-1)],TairsMeans[(t+1)])
  #     }
  #   }
  # }
  #
  # #TairsVar[is.infinite(TairsVar)] <- 100
  # TairsPrec <- 1/TairsVar
  # TairsPrec[is.infinite(TairsPrec)] <- 1/100

  #dat <- list(Sf=SfsMeans,Sfprec=1/SfsVar)
  dat <- list(TairMu=TairsMeans,TairPrec=TairsPrec)
  #print("finished createTairs")
  return(dat)
}
