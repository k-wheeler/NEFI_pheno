##' Creates data object with the desired Sf means and variances
##'
##' @param latitude The desired latitude if not willowCreek
##' @param longitude The desired longitude if not willowCreek
##' @param siteName The site name
##' @param dataDirectory The data directory
##' @param endDate The end date
##' @param GEFS_Files A vector of the file names for the GEFS files
##' @param GEFS_Directory The directory where the GEFS files are located
##' @export
createSf <- function(lat="",long="",years,siteName,dataDirectory,endDate,GEFS_Files,GEFS_Directory) {
  calDates <- seq(as.Date(paste(years[1],"-01-01",sep="")),as.Date(paste((years[length(years)]-1),"-12-31",sep="")),"day")
  #calDates <- calDates[as.numeric(format(calDates,"%j"))%in% seq(1,181)]

  ##The sources of the calibration and current measurements differ between willowCreek and other sites
  if(siteName=="willowCreek"){
    TairsCalInd <- download_US_WCr_met(start_date=calDates[1],end_date=calDates[length(calDates)])
    ##Need to pad TairsCal to have ensembles
    TairsCal <- matrix(ncol=10,nrow=length(TairsCalInd))
    for(e in 1:10){
      TairsCal[,e] <- TairsCalInd
    }
    TairsCurrentInd <- download_US_WCr_met(start_date=as.Date("2019-01-01"),end_date=endDate)
  }
  else{
    TairsCal <- load_ERA5_Tair(lat=lat,long=long,years=years,dataDirectory=dataDirectory) ##columns are each an ensemble (not divided by year)
    TairsCurrentInd <- load_NOAA_met(lat=lat,long=long,years=years) ##Array of numeric values
  }

  ##GEFS Forecast (same for all sites) and pad TairsCurrent to be ensembles
  TairsCurrent <- matrix(ncol=length(GEFS_Files),nrow=length(TairsCurrentInd))
  TairsForecast <- matrix(ncol=length(GEFS_Files),nrow=15)
  for(e in 1:length(GEFS_Files)){
    TairsForecastInd <- load_GEFS_Forecast(dataDirectory=GEFS_Directory,fileName=GEFS_Files[e])
    TairsForecast[,e] <- TairsForecastInd
    TairsCurrent[,e] <- TairsCurrentInd
  }

  ##Create Sfs
  SfsALL <- matrix(nrow=0,ncol=length(calDates))
  for(e in 1:ncol(TairsCal)){
    Sfs <- calSf(Tairs=TairsCal[,e],dates=calDates)
    SfsALL <- rbind(SfsALL,Sfs)
  }
  ##Current year
  curDates=seq(as.Date("2019-01-01"),endDate,"day")
  for(e in 1:length(GEFS_Files)){
    curDates <- na.omit(curDates)
    Tairs <- c(TairsCurrent[,e],TairsForecast[,e])[1:length(curDates)]
    Sfs <- calSf(Tairs=Tairs,dates=curDates)
    Sfs <- c(Sfs,rep(NA,(365-length(Sfs)))) ##Done to make all years 365 days
    SfsALL <- rbind(SfsALL,Sfs)
  }

    SfsMeans <- colMeans(SfsALL)
    SfsVar <- apply(SfsALL,MARGIN=2,FUN=var)
    SfsVar[SfsVar==0] <- 0.001
    dat <- list(Sf=SfsMeans,Sfprec=1/SfsVar)
    return(dat)
}
