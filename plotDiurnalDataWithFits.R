#!/usr/bin/env Rscript

#install.packages("devtools")
#library("devtools")
#install_github("EcoForecast/ecoforecastR")
library("ecoforecastR")
library("rjags")

diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  #print(k)
  bk <- which(round(xseq,digits=1)==k)
  #print(bk)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  #print(length(c(left,right)))
  return(c(left,right))
}

siteName <- "harvard"
#diurnalFiles <- intersect(dir(path="dailyNDVI_GOES",pattern="varBurn2.RData"),dir(path="dailyNDVI_GOES",pattern=siteName))
#diurnalFiles <- dir(path="dailyNDVI_GOES",pattern=siteName)
diurnalFiles <- intersect(dir(path="PhenologyForecastData/",pattern=siteName),dir(path="PhenologyForecastData/",pattern="GOES_NDVI_Diurnal"))
xseq <- seq(0,25,0.1)
#i=1
outputFileName <- paste(siteName,"_DiurnalFits_withData.pdf",sep="")
pdf(file=outputFileName,width=45,height=40)
par(mfrow=c(5,5))
for(i in 1:length(diurnalFiles)){
  #dayData <- read.csv(paste("dailyNDVI_GOES/",diurnalFiles[i],sep=""),header=FALSE)
  dayData <- read.csv(paste("PhenologyForecastData/",diurnalFiles[i],sep=""),header=FALSE)
  yearDay <- substr((strsplit(diurnalFiles[i],"_")[[1]][4]),1,7)
  print(diurnalFiles[i])

  plot(as.numeric(dayData[3,]),as.numeric(dayData[2,]),main=diurnalFiles[i],ylim=c(0,1),xlim=c(0,25))
  #fitFileName <- paste(siteName,"_",day,"_varBurn2.RData",sep="")
  fitFileName <- paste("PhenologyForecastData/GOES_DiurnalFits/",siteName,"_",yearDay ,"_varBurn.RData",sep="")
  if(file.exists(fitFileName)){
    load(fitFileName)
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
    abline(v=12,col="red")
  }
}
dev.off()
