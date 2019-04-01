##' Creates data object with the desired Sf means and variances
##'
##' @param lat The desired latitude if not willowCreek
##' @param long The desired longitude if not willowCreek
##' @param dates The desired dates
##' @param siteName The site name
##' @param dataDirectory The data directory
##' @param endDate The end date
##' @param GEFS_Files A vector of the file names for the GEFS files
##' @param GEFS_Directory The directory where the GEFS files are located
##' @export
createSf <- function(lat="",long="",dates,siteName,dataDirectory,endDate,GEFS_Files,GEFS_Directory,forecastLength) {
  years <- lubridate::year(dates)
  calDates <- seq(dates[1],as.Date(paste((years[length(years)]-1),"-12-31",sep="")),"day")

  ##The sources of the calibration and current measurements differ between willowCreek and other sites
  if(siteName=="willowCreek"){
    TairsCalInd <- download_US_WCr_met(start_date=calDates[1],end_date=calDates[length(calDates)])

    ##Need to pad TairsCal to have ensembles
    TairsCal <- matrix(ncol=10,nrow=length(TairsCalInd))
    for(e in 1:10){
      TairsCal[,e] <- TairsCalInd
    }
    TairsCurrentInd <- download_US_WCr_met(start_date=as.Date("2019-01-01"),end_date=(endDate-forecastLength))

  }
  else{
    TairsCal <- load_ERA5_Tair(lat=lat,long=long,years=seq(years[1],2018)) ##columns are each an ensemble (not divided by year)
    TairsCurrentInd <- load_NOAA_met(lat=lat,long=long,years=seq(2019,2019)) ##Array of numeric values
  }

  ##GEFS Forecast (same for all sites) and pad TairsCurrent to be ensembles
  TairsCurrent <- matrix(ncol=length(GEFS_Files),nrow=length(TairsCurrentInd))
  TairsForecast <- numeric()
  for(e in 1:length(GEFS_Files)){
    TairsForecastInd <- load_GEFS_Forecast(dataDirectory=GEFS_Directory,fileName=GEFS_Files[e])
    TairsForecast <- cbind(TairsForecast,TairsForecastInd)
    TairsCurrent[,e] <- TairsCurrentInd
  }

  ##Create Sfs
  SfsALL <- matrix(nrow=0,ncol=length(calDates))

  for(e in 1:ncol(TairsCal)){
    Sfs <- calSf(Tairs=TairsCal[,e],dates=calDates)
    SfsALL <- rbind(SfsALL,Sfs)
  }
  SfsMeansCal <- colMeans(SfsALL)
  SfsVarCal <- apply(SfsALL,MARGIN=2,FUN=var)
  SfsVarCal[SfsVarCal==0] <- 0.001

  ##Current year
  curDates=seq(as.Date("2019-01-01"),endDate,"day") ##Includes forecasted period

  SfsALL <- matrix(nrow=0,ncol=length(curDates))
  for(e in 1:length(GEFS_Files)){
    Tairs <- c(TairsCurrent[,e],TairsForecast[,e])
    print("Tairs:")
    print(Tairs)
    Sfs <- calSf(Tairs=Tairs,dates=curDates)
    print("Sfs:")
    print(Sfs)
    SfsALL <- rbind(SfsALL,Sfs)
  }

  SfsMeansCur <- colMeans(SfsALL)
  SfsVarCur <- apply(SfsALL,MARGIN=2,FUN=var)
  SfsVarCur[SfsVarCur==0] <- 0.001
  SfsMeans <- c(SfsMeansCal,SfsMeansCur)
  SfsVar <- c(SfsVarCal,SfsVarCur)

  dat <- list(Sf=SfsMeans,Sfprec=1/SfsVar)
  return(dat)
}
