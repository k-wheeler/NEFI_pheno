#' Main function to calculate all NDVI values for a day for a list of sites
#'
#' @param day The desired day of year
#' @param siteData matrix of site data where the sites are in individual rows and the columns are siteName, latitude, longitude, timezone (Note: all timezone values must be the same)
#' @param year The desired year
#' @param TZ The timezone of the sites in siteData (e.g. eastern US sites have a TZ value of 5).
#' @param dataPath The directory where all of the GOES data is located
#' @param TZ_name The name of the time zone (e.g., "America/New_York")
#' @param savePath The directory where you want to save the output
#' @param mVersion The mode/version of the ABI (3,4, or 6)
#' @import suncalc
#' @import lubridate
#' @export
calculateNDVI_GOES_MAIN <- function(day,siteData,year,TZ,dataPath,TZ_name,savePath,mVersion){
  date.val <- as.Date(as.numeric(day),origin=as.Date(paste(as.character(as.numeric(year)-1),"-12-31",sep="")))
  print(date.val)
  if(year==2017 && day < 321){
    orbitVersion <- "OLD"
  }
  else{
    orbitVersion <- "NEW"
  }

  filestrACM <- paste("OR_ABI-L2-ACMC-M",mVersion,"_G16_s",year,day,sep="")
  ACM.files <- dir(path=dataPath,pattern=filestrACM)

  #if(!dir.exists((paste(dataPath,dir(path=dataPath,pattern=filestrACM),sep="")))){
  if(length(ACM.files>1)){
    print("ACM.files>1")
    day.time.vals <- character()
    times <- character()
    NDVI.vals <- matrix(ncol=nrow(siteData),nrow=0)
    for(j in 1:length(ACM.files)){
      print(j)
      day.time <- substr(ACM.files[j],24,34)
      day.time.vals <- c(day.time.vals,day.time)
      hr <- as.numeric(substr(day.time,8,9))-TZ
      mt <- substr(day.time,10,11)
      times <- c(times,(as.numeric(hr)+as.numeric(mt)/60))
      print(paste("Time:",(as.numeric(hr)+as.numeric(mt)/60)))
      newNDVI <- createNDVI_sub(siteData=siteData,orbitVersion=orbitVersion,day.time=day.time,dataPath=dataPath,mVersion=mVersion)
      NDVI.vals <- rbind(NDVI.vals,newNDVI)
    }
    for(i in 1:nrow(siteData)){
      siteName <- as.character(siteData[i,1])
      fileName <- paste(savePath,"GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="")
      output <- rbind(day.time.vals,NDVI.vals[,i],times)
      suntimes <- getSunlightTimes(date=date.val,lat=as.numeric(siteData[i,2]),lon=as.numeric(siteData[i,3]),keep=c("nauticalDawn","nauticalDusk"),tz = TZ_name)
      dawnTime <- lubridate::hour(suntimes$nauticalDawn)+(lubridate::minute(suntimes$nauticalDawn)/60)
      duskTime <- lubridate::hour(suntimes$nauticalDusk)+(lubridate::minute(suntimes$nauticalDusk)/60)

      output[2,as.numeric(output[3,])<as.numeric(dawnTime)] <- NA
      output[2,as.numeric(output[3,])>as.numeric(duskTime)] <- NA
      write.table(output,file=fileName,sep=",",col.names=FALSE,row.names=FALSE)
    }
  }
  else{
    print("No ACM files")
    #createEmptyFiles(siteData=siteData,day=day,year=year)
  }
}
# else{
#   print("File already exists")
#   #createEmptyFiles(siteData=siteData,day=day,year=year)
# }
# }
