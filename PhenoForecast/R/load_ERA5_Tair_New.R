##' Calculates cummulative Tair (within one year for one ensemble)
##'
##' @param lat The site latitude
##' @param long The site longitude
##' @param endDate
##' @param calDatesT
##' @param ERA5dataFolder
##' @param TZ_name
##' @import xts
##' @import ncdf4
##' @export
load_ERA5_Tair_New <- function(lat="",long="",endDate="",calDatesT=TRUE,ERA5dataFolder,TZ_name="America/New_York") {
  if(calDatesT){ #Need to include calibration data
    calFileName <- dir(path=ERA5dataFolder,pattern=paste(endDate,"_era5AirTemperatureMembers",sep=""))

    if(length(calFileName)==0){
      downloadERA5Calibration(var="ensemble_members") ##***Need to add this function***
      calFileName <- dir(path=ERA5dataFolder,pattern=paste(endDate,"_era5AirTemperatureMembers",sep=""))
    }
    #print(calFileName)
    #Load data
    ensembleFile <- nc_open(paste(ERA5dataFolder,calFileName,sep=""))

    Tairs <- ncvar_get(ensembleFile)-273 #Convert from Kelvin to C


    timeHours <- ensembleFile$dim$time$vals #Hours since 1900-01-01 00:00:00.0

    ##Convert times to actual times
    times <- as.POSIXct(timeHours*3600, origin = "1900-01-01",tz = "GMT")
    attributes(times)$tzone <- TZ_name
    #Daily average
    allDates <- lubridate::date(times)
    dates <- seq(lubridate::date(times[5]),lubridate::date(times[length(times)]),"day")
    TairsDaily <- matrix(nrow=10,ncol=length(dates))
    print(length(dates))
    print(dim(Tairs))
    #print(colnames(Tairs))
    print(head(Tairs[,1,]))
    print(head(Tairs[,2,]))
    for(d in 1:length(dates)){
      subTairs <- Tairs[,allDates==dates[d]]
      TairsDaily[,d] <- apply(subTairs,MARGIN=1,mean)
    }
  }
  ##Current Year (already downloaded)
    ##Will fill in once I get the calibration done
  return(TairsDaily)

}
library(ncdf4)
ERA5dataFolder="/projectnb/dietzelab/kiwheel/ERA5/Data/dukehw/"
TZ_name="America/New_York"

lat=""
long=""
siteName="dukehw"
endDate=""


