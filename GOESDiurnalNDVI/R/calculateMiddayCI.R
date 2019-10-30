#' Calculate credible interval on midday NDVI value
#'
#' @param siteName The site name used for file naming
#' @param year The desired year
#' @param day The desired day of year
#' @param savePath The folder where the model output is saved
#' @export
#' @import rjags
#' @import runjags
calculateMiddayCI <- function(siteName,year,day,savePath){
  fileName <- paste(savePath,siteName,"_",year,day,"_varBurn.RData",sep="")
  load(fileName)
  dat <- read.csv(paste(savePath,"GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep=""),header=FALSE)
  #print(dim(dat))
  xseq <- seq(dat[3,1],dat[3,ncol(dat)],0.001)
  out.mat <- as.matrix(var.burn)
  a <- out.mat[,1]
  c <- out.mat[,2]
  k <- out.mat[,3]
  ycred <- matrix(0,nrow=10000,ncol=length(xseq))
  for(g in 1:10000){
    Ey <- diurnalExp(a=a[g],c=c[g],k=k[g],xseq=xseq)
    ycred[g,] <- Ey
  }
  ci <- apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE)
  return(ci)
}

