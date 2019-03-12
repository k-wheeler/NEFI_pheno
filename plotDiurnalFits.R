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

outputFileName <- "ALL_DiurnalFits.pdf"
#pdf(file=outputFileName,width=45,height=40)
xseq <- seq(0,25,0.1)
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
#print(siteData[,1])
iseq <- c(13,14)
for(s in iseq){
  print(s)
  siteName <- as.character(siteData[s,1])
  print(siteName)
  outputFileName <- paste(siteName,"_ALL_DiurnalFits3.pdf",sep="")
  pdf(file=outputFileName,width=45,height=40)
  par(mfrow=c(5,5))
  diurnalFiles <- intersect(dir(pattern="varBurn4.RData"),dir(pattern=siteName))
  
  for(i in 1:length(diurnalFiles)){
    load(diurnalFiles[i])
    dy <- strsplit(diurnalFiles[i],"_")[[1]][2]
    print(diurnalFiles[i])
    print(dy)
    print(as.numeric(dy))
    if(as.numeric(dy)<182){
      yr <- "2018"
    }
    else{
      yr <- "2017"
    }
    dayDataFile <- paste("dailyNDVI_GOES/",siteName,"_GOES_diurnal_",yr,dy,".csv",sep="")
    if(file.exists(dayDataFile)){
      dat <- read.csv(dayDataFile,header=FALSE)
    }
    else{
      dat <- matrix(ncol=5,nrow=3)
      dat[1,] <- c(5,6,7,8,9)
      dat[2,] <- c(NA,NA,NA,NA,NA)
      dat[3,] <- c(5,6,7,8,9)
    }
    if(typeof(var.burn)==typeof(FALSE)){
      print(paste(diurnalFiles[i], " did not converge",sep=""))
      #plot(as.numeric(dat[3,]),as.numeric(dat[2,]),main=paste("Didn't",diurnalFiles[i],sep=" "),xlab="Time",ylab="NDVI",ylim=c(0,1),xlim=c(0,25))
    }
    else{
      out.mat <- as.matrix(var.burn)
      a <- out.mat[,1]
      rndNums <- sample(1:length(a),10000,replace=T)
      a <- a[rndNums]
      c <- out.mat[rndNums,2]
      k <- out.mat[rndNums,3]
      ycred <- matrix(0,nrow=10000,ncol=length(xseq))
      for(g in 1:10000){
        Ey <- diurnalExp(a=a[g],c=c[g],k=k[g],xseq=xseq)
        ycred[g,] <- Ey
      }
      ci <- apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE)
      #plot(x=list(),y=list(),main=diurnalFiles[i],xlab="Time",ylab="NDVI",ylim=c(0,1),xlim=c(0,25))
      if(length(na.omit(as.numeric(dat[2,])))>0){
        plot(as.numeric(dat[3,]),as.numeric(dat[2,]),main=diurnalFiles[i],xlim=c(0,25),pch=20,cex=2,ylim=c(0,1))
        ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
        lines(xseq,ci[2,],col="black")
        points(as.numeric(dat[3,]),as.numeric(dat[2,]),pch=20,cex=2)
        abline(v=12,col="red")
      }
    }
  }
  dev.off()
}
