#' Plots the fitted curve with credible interval for a day
#'
#' @param siteName The site name used for file naming
#' @param year The desired year
#' @param day The desired day of year
#' @param savePath The folder where the model output is saved
#' @export
#' @import ecoforecastR
plotCI <- function(siteName,year,day,savePath){
  ci <- calculateMiddayCI(siteName=siteName,year=year,day=day,savePath=savePath)
  date <- as.Date(as.numeric(day),origin=as.Date(paste(as.character(as.numeric(year)-1),"-12-31",sep="")))
  dataFileName <- paste(savePath, "GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="")
  dat <- read.csv(dataFileName,header=FALSE)
  plot(x=list(),y=list(),main=date,ylim=c(0,1),xlim=c(0,24),ylab="NDVI",xlab="Hour",cex=2.5)
  #polygon(x=c(10,14,14,10),y=c(-1,-1,1.2,1.2),col="lightgray",border=NA)
  xseq <- seq(dat[3,1],dat[3,ncol(dat)],0.001)
  ecoforecastR::ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
  lines(xseq,ci[2,],col="black")

  points(as.numeric(dat[3,]),as.numeric(dat[2,]),pch=20)
  #abline(v=12,col="red")
}
