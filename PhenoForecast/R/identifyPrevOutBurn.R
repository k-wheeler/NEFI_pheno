#' Identifies the most recent previous forecast file
#'
#' @param siteName
#' @param forecastStr
#' @param startDate
#' @param endDate
#' @param dataDirectory
#'
#' @return
#' @export
#'
#' @examples
identifyPrevOutBurn <- function(siteName,forecastStr,startDate,endDate,dataDirectory){
  potentialDates <- seq(as.Date("2019-10-01"),(endDate-1),"day")
  fileFound <- FALSE
  i <- length(potentialDates)
  while(!fileFound){
    fileName <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",potentialDates[i],"/",
                      siteName,"_",startDate,"_",endDate,"_",forecastStr,
                      "_outBurn.RData",sep="")
    print(fileName)
    if(file.exists(fileName)){
      fileFound <- TRUE
    }else if(i==1){
      fileFound <- TRUE
      fileName <- NA
    }else{
      i <- i - 1
    }
  }
  print("fileName:")
  print(fileName)
  return(fileName)
}
