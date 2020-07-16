##' For PhenoCam data, construct the data object for input into MCMC
##'
##' @param URL PhenoCam network URL
##' @param siteName The site name
##' @param startDate The start date
##' @param endDate The end date
##' @param metric The desired metric ("GCC" or "RCC")
##' @export
PC_data <- function(siteName,URL,startDate,endDate,seasonOrder="AS",metric="GCC") {
  ##Data
  if(metric=="GCC"){
    fileName <- paste(siteName,"_",startDate,"_",endDate,"PC.RData",sep="")
  }else if(metric=="RCC"){
    fileName <- paste(siteName,"_",startDate,"_",endDate,"PC_RCC.RData",sep="")
  }

  years <- seq(lubridate::year(startDate),lubridate::year(endDate))
  if(!file.exists(fileName)){
    ##Download all data for those years
    PC.data <- subset(download.phenocam(URL),year%in%years)
    ##Index for the startDate and endDate
    PC.startDayIndex <- which(as.Date(PC.data$date)==startDate)
    PC.endDayIndex <- which(as.Date(PC.data$date)==endDate)
    ##Subset data for specific date range
    PC.data <- PC.data[PC.startDayIndex:PC.endDayIndex,]
    PC.time <-  as.Date(PC.data$date)
    if(metric=="GCC"){
      y <- PC.data$gcc_mean
      obs.prec <- 1/(PC.data$gcc_std**2)
    }else if(metric=="RCC"){
      y <- PC.data$rcc_mean
      obs.prec <- 1/(PC.data$rcc_std**2)
    }

    #print(PC.data$gcc_std)

    obs.prec[is.infinite(obs.prec)] <- 1/(0.001**2)
    obs.prec[is.na(obs.prec)] <- 1/(0.001**2)
    x <- lubridate::yday(PC.time)
    ##If the season order is autumn and then spring, it adds 365 to DOY of spring
    if(seasonOrder=="AS"){
      bk <- which(x==365)+1
      for(i in bk:length(x)){
        x[i] <- x[i]+365
      }
    }
    #print(obs.prec)
    #print(length(x))
    #print(length(obs.prec))
    data <- list(x=x,y=y,n=length(y),obs.prec=obs.prec)
    save(data,file=fileName)
  }
  load(fileName)
  return(data)
}

