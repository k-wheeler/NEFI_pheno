##' Function to plot phenology forecast output
##'
##' @param siteName The site name to be printed on the graph
##' @param forecastType The type of forecast (randomWalk or logistic)
##' @param URL The PhenoCam URL
##' @param forecastLength The number of days in the future you want to forecast
##' @param out.mat The predict variables of the MCMC output of the forecast in matrix form
##' @export
##' @import ecoforecastR
##' @import rjags
plotForecastOutput <- function(siteName,forecastType,URL,forecastLength,out.mat,days){
  ##Download the phenocam data
  phenoData <- matrix(nrow=0,ncol=32)
  print(URL)
  for(u in 1:length(URL)){
    print(URL[u])
    phenoDataSub <- download.phenocam(URL[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }
  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  phenoData <- phenoData[phenoData$date<endDate,]

  p <- phenoData$gcc_mean
  time <-  as.Date(phenoData$date)

  ###Pad the x (time) and y (GCC) with future days and GCC values of NA to stimulate a forecast
  p <- c(p,rep(x=NA,times=forecastLength)) #Padded with NA's to forecast for one month into the future
  timeForecast <- c(time,seq.Date(from=time[length(time)],by="day",length.out=forecastLength))

  #p <- rescaleObs(time=timeForecast,vals=p) ##Rescale phenocam data between (0,1)

  ci <- apply(out.mat,2,quantile,c(0.025,0.5,0.975)) #Computes the 95% credible interval (CI)
  ##Plot
  print(length(ci[2,]))
  print(length(days))
  plot(days,ci[2,],type='n',xlab="Time",ylab="Percent Canopy",main=paste(siteName,forecastType),cex.lab=1.5,cex.main=2,ylim=c(0,1))

  ciEnvelope(days,ci[1,],ci[3,],col="lightBlue")
  #lines(days,ci[1,],col="purple")
  #lines(days,ci[3,],col="green")
  lines(days,ci[2,],col="black")
  #points(timeForecast,p,pch="+",cex=0.5)
  abline(v=time[length(time)],col="red")

}
