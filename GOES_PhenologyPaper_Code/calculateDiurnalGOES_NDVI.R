#' calculateDiurnalGOES_NDVI
#'
#' @param savePath
#' @param allSiteData
#' @param GOESdataPath
#' @param calculatedNDVIGOESpath
#' @param year
#' @param TZ
#' @param TZ_name
#' @return
#' @export
#' @import GOESDiurnalNDVI
#'
#' @examples
calculateDiurnalGOES_NDVI <- function(savePath,allSiteData,GOESdataPath,calculatedNDVIGOESpath,year,TZ,TZ_name){
  # library("doParallel")
  # n.cores <- 16 #Included if you want and are able to run in Parallel using the doParallel package 
  # registerDoParallel(cores=n.cores)
  days <- seq(1,365)
  if(year==2018){
    mVersions <- rep(3,365) #Specific for 2018
  }else if(year==2019){
    mVersions <- c(rep(3,92),rep(6,273)) #Specific for 2019
  }else{
    print("Year not recognized")
  }

  siteData <- allSiteData[allSiteData$TZ==TZ,]
  for(i in 1:length(days)){
    # output <- #Included if you want to run in parallel with the doParallel package
    # foreach(i=1:length(days)) %dopar% { 
    mVersion <- mVersions[i]
    day <- days[i]
    print(day)
    calculateNDVI_GOES_MAIN(day=day,siteData=siteData,year=year,TZ=TZ,
                            dataPath=GOESdataPath,TZ_name=TZ_name,
                            savePath=calculatedNDVIGOESpath,mVersion = mVersion)
  }
}
