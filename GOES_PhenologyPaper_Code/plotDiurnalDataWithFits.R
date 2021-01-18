#' diurnalExp
#'
#' @param a 
#' @param c 
#' @param k 
#' @param xseq 
#'
#' @return
#' @export
#'
#' @examples
diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  bk <- which(round(xseq,digits=1)==k)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  return(c(left,right))
}

#' plotDiurnalDataWithFits
#'
#' @param siteName Desired site name based (e.g., "harvard")
#' @param yr Desired year (e.g., 2019)
#' @return
#' @import ecoforecastR
#' @import rjags
#' @export
#'
#' @examples
plotDiurnalDataWithFits <- function(calculatedNDVIGOESpath,DiurnalFitSavePath,siteName,yr){
  outputFileName <- paste(siteName,"_DiurnalFits_withData.pdf",sep="")
  diurnalFiles <- intersect(dir(path=paste(calculatedNDVIGOESpath,"GOES_NDVI_DiurnalData",sep=""),pattern=siteName),dir(path=paste(savePath,"GOES_NDVI_DiurnalData",sep=""),pattern=as.character(yr)))
  
  xseq <- seq(0,25,0.1)

  pdf(file=outputFileName,width=25,height=25)
  par(mfrow=c(5,5))
  for(i in seq(1,length(diurnalFiles),5)){
    dayData <- read.csv(paste(calculatedNDVIGOESpath,diurnalFiles[i],sep=""),header=FALSE)
    yearDay <- substr((strsplit(diurnalFiles[i],"_")[[1]][4]),1,7)
    print(diurnalFiles[i])
    
    plot(as.numeric(dayData[3,]),as.numeric(dayData[2,]),main=diurnalFiles[i],ylim=c(0,1),xlim=c(0,25))

    fitFileName <- intersect(dir(path=DiurnalFitSavePath,pattern=siteName),dir(path=savePath,pattern=yearDay))
    print(fitFileName)
    if(length(fitFileName)>0){
      load(paste(DiurnalFitSavePath,fitFileName,sep=""))
      out.mat <- as.matrix(var.burn)
      rndNums <- sample.int(nrow(out.mat),10000,replace=T)
      a <- out.mat[rndNums,1]
      c <- out.mat[rndNums,2]
      k <- out.mat[rndNums,3]
      ycred <- matrix(0,nrow=10000,ncol=length(xseq))
      for(g in 1:10000){
        Ey <- diurnalExp(a=a[g],c=c[g],k=k[g],xseq=xseq)
        ycred[g,] <- Ey
      }
      ci <- apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE)
      ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
      lines(xseq,ci[2,],col="black")
      points(as.numeric(dayData[3,]),as.numeric(dayData[2,]))
    }
  }
  dev.off()
}

