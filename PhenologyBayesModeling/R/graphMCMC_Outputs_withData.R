##' Plot the MCMC outputs with data (rough and not generalized)
##'
##' @param outputFileName The filename for the graphs
##' @param siteFileName The file name of the site data
##' @param iseq The row numbers for the desired sites in the siteFileName
##' @param startDaySH The start day for the shrubland sites
##' @param endDaySH The end day for the shrubland sites
##' @param startDayDB The start day for the deciduous broadleaf sites
##' @param endDayDB The end day for the deciduous broadleaf sites
##' @import rjags
##' @export
graphMCMC_Outputs_withData <- function(outputFileName,siteFileName,iseq,startDaySH,endDaySH,startDayDB,endDayDB){
  pdf(file=outputFileName,width=8,height=8)
  siteData <- read.csv(siteFileName,header=TRUE)

  for (i in iseq){
    siteName <- as.character(siteData$siteName[i])
    print(siteName)
    URL <- as.character(siteData$URL[i])
    PFT <- as.character(siteData$PFT[i])
    lat <- as.character(siteData$Lat[i])
    long <- as.character(siteData$Long[i])
    if(PFT=="DB"){
      startDay <- startDayDB
      endDay <- endDayDB
    }
    else if(PFT=="SH"){
      startDay <- startDaySH
      endDay <- endDaySH
    }
    xseq <- seq(startDay,endDay,1)
    #GOES
     data.GOES = GOES_data(siteName=siteName,startDay = startDay,endDay = endDay,lat=lat,long=long)
     inputFileName <- paste(siteName,"_GOES_varBurn.RData",sep="")
     load(inputFileName)

     var.mat<-as.matrix(GOES.md.out)
     GOES.c <- var.mat[,3]
     GOES.d <- var.mat[,4]
     if(PFT=="DB"){
       GOES.Tran1<-var.mat[,1]
       GOES.Tran2<-var.mat[,2]
       GOES.c <- var.mat[,5]
       GOES.d <- var.mat[,6]

     }
     else if(PFT=="SH"){
       GOES.Tran1<-var.mat[,1]
       GOES.Tran2<-var.mat[,5]
       GOES.c <- var.mat[,3]
       GOES.d <- var.mat[,4]
     }
     ci.GOES <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)
    # print("GOES Done")

    #MODIS NDVI:
    data.MODIS.N = MODIS_data(siteName=siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="NDVI")
    inputFileName <- paste(siteName,"_MODIS_NDVI_varBurn.RData",sep="")
    load(inputFileName)

    var.mat<-as.matrix(MODIS.N.md.out)

    if(PFT=="DB"){
      MODIS.N.Tran1<-var.mat[,1]
      MODIS.N.Tran2<-var.mat[,2]
      MODIS.N.c <- var.mat[,5]
      MODIS.N.d <- var.mat[,6]

    }
    else if(PFT=="SH"){
      MODIS.N.Tran1<-var.mat[,1]
      MODIS.N.Tran2<-var.mat[,5]
      MODIS.N.c <- var.mat[,3]
      MODIS.N.d <- var.mat[,4]
    }

    ci.MODIS.N <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)
    print("MODIS NDVI Done")

    #MODIS EVI:
    data.MODIS.E = MODIS_data(siteName=siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="EVI")
    inputFileName <- paste(siteName,"_MODIS_EVI_varBurn.RData",sep="")
    load(inputFileName)

    var.mat<-as.matrix(MODIS.E.md.out)

    if(PFT=="DB"){
      MODIS.E.Tran1<-var.mat[,1]
      MODIS.E.Tran2<-var.mat[,2]
      MODIS.E.c <- var.mat[,5]
      MODIS.E.d <- var.mat[,6]

    }
    else if(PFT=="SH"){
      MODIS.E.Tran1<-var.mat[,1]
      MODIS.E.Tran2<-var.mat[,5]
      MODIS.E.c <- var.mat[,3]
      MODIS.E.d <- var.mat[,4]
    }

    ci.MODIS.E <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)
    print("MODIS EVI Done")

    # #PC:
     data.PC = PC_data(siteName=siteName,URL=URL,startDay=startDay,endDay=endDay)
     inputFileName <- paste(siteName,"_PC_varBurn.RData",sep="")
     load(inputFileName)

     var.mat<-as.matrix(PC.md.out)
     if(PFT=="DB"){
       PC.Tran1<-var.mat[,1]
       PC.Tran2<-var.mat[,2]
       PC.c <- var.mat[,5]
       PC.d <- var.mat[,6]

     }
     else if(PFT=="SH"){
       PC.Tran1<-var.mat[,1]
       PC.Tran2<-var.mat[,5]
       PC.c <- var.mat[,3]
       PC.d <- var.mat[,4]
     }
     ci.PC <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)
     print("PC Done")

    par(mfrow=c(1,1))
    plot(x=list(),y=list(),xlim=c(100,550),ylim=c(-0.2,1.2),ylab="Value",xlab="Day of Year",main=paste(siteName,"GOES"),cex.axis=2,cex.lab=2,cex.main=2)
    lines(xseq,ci.GOES[2,],col="black",lwd=2)
    lines(xseq,ci.GOES[1,],col="black", lty = 2,lwd=2)
    lines(xseq,ci.GOES[3,],col="black", lty = 2,lwd=2)
    abline(v=mean(GOES.Tran1),col="black")
    abline(v=mean(GOES.Tran2),col="black")
    points(data.GOES$x,rescale(c=mean(GOES.c),d=mean(GOES.d),yseq=data.GOES$y),col="Black",pch=20)

    plot(x=list(),y=list(),xlim=c(100,550),ylim=c(-0.2,1.2),ylab="Value",xlab="Day of Year",main=paste(siteName,"PhenoCam"),cex.axis=2,cex.lab=2,cex.main=2)
    lines(xseq,ci.PC[2,],col="cyan",lwd=2)
    lines(xseq,ci.PC[1,],col="cyan", lty = 2,lwd=2)
    lines(xseq,ci.PC[3,],col="cyan", lty = 2,lwd=2)
    abline(v=mean(PC.Tran1),col="cyan")
    abline(v=mean(PC.Tran2),col="cyan")
    points(data.PC$x,rescale(c=mean(PC.c),d=mean(PC.d),yseq=data.PC$y),col="cyan",pch=20)

    plot(x=list(),y=list(),xlim=c(100,550),ylim=c(-0.2,1.2),ylab="Value",xlab="Day of Year",main=paste(siteName,"MODIS NDVI"),cex.axis=2,cex.lab=2,cex.main=2)
    lines(xseq,ci.MODIS.N[2,],col="red",lwd=2)
    lines(xseq,ci.MODIS.N[1,],col="red", lty = 2,lwd=2)
    lines(xseq,ci.MODIS.N[3,],col="red", lty = 2,lwd=2)
    abline(v=mean(MODIS.N.Tran1),col="red")
    abline(v=mean(MODIS.N.Tran2),col="red")
    points(data.MODIS.N$x,rescale(c=mean(MODIS.N.c),d=mean(MODIS.N.d),yseq=data.MODIS.N$y),col="red",pch=20)

    plot(x=list(),y=list(),xlim=c(100,550),ylim=c(-0.2,1.2),ylab="Value",xlab="Day of Year",main=paste(siteName,"MODIS EVI"),cex.axis=2,cex.lab=2,cex.main=2)
    lines(xseq,ci.MODIS.E[2,],col="red",lwd=2)
    lines(xseq,ci.MODIS.E[1,],col="red", lty = 2,lwd=2)
    lines(xseq,ci.MODIS.E[3,],col="red", lty = 2,lwd=2)
    abline(v=mean(MODIS.E.Tran1),col="red")
    abline(v=mean(MODIS.E.Tran2),col="red")
    points(data.MODIS.E$x,rescale(c=mean(MODIS.E.c),d=mean(MODIS.E.d),yseq=data.MODIS.E$y),col="red",pch=20)

  }
  dev.off()
}
