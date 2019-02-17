##' Plot the MCMC outputs (rough and not generalized)
##'
##' @param outputFileName The filename for the graphs
##' @param siteFileName The file name of the site data
##' @param iseq The row numbers for the desired sites in the siteFileName
##' @param startDay The start day
##' @param endDay The end day
##' @import rjags
##' @export
graphMCMC_Outputs <- function(outputFileName,siteFileName,iseq,startDay,endDay){
  pdf(file=outputFileName,width=8,height=8)
  siteData <- read.csv(siteFileName,header=FALSE)
  xseq <- seq(startDay,endDay,1)

  for (i in iseq){
    siteName <- as.character(siteData[i,1])
    print(siteName)
    URL <- as.character(siteData[i,4])
    PFT <- as.character(siteData[i,5])

    #GOES
    #data.GOES = GOES_data(siteName=siteName,startDay = startDay,endDay = endDay)
    inputFileName <- paste(siteName,"_GOES_varBurnINITS.RData",sep="")
    load(inputFileName)

    var.mat<-as.matrix(GOES.md.out)
    if(PFT=="DB"){
      GOES.a <- var.mat[,1]
      GOES.b <- var.mat[,2]
      GOES.trans <- mean(-(GOES.a/GOES.b))
    }
    else if(PFT=="SH"){
      GOES.trans <- mean(var.mat[,5])
    }
    ci.GOES <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)

    #MODIS:
    #data.MODIS = MODIS_data(siteName=siteName)
    inputFileName <- paste(siteName,"_MODIS_varBurnINITS.RData",sep="")
    load(inputFileName)

    var.mat<-as.matrix(MODIS.md.out)
    if(PFT=="DB"){
      MODIS.a <- var.mat[,1]
      MODIS.b <- var.mat[,2]
      MODIS.trans <- mean(-(MODIS.a/MODIS.b))
    }
    else if(PFT=="SH"){
      MODIS.trans <- mean(var.mat[,5])
    }
    ci.MODIS <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)

    #PC:
    #data.PC = PC_data(URL=URL,startDay=startDay,endDay=endDay)
    inputFileName <- paste(siteName,"_PC_varBurnINITS.RData",sep="")
    load(inputFileName)

    var.mat<-as.matrix(PC.md.out)
    if(PFT=="DB"){
      PC.a <- var.mat[,1]
      PC.b <- var.mat[,2]
      PC.trans <- mean(-(PC.a/PC.b))
    }
    else if(PFT=="SH"){
      PC.trans <- mean(var.mat[,5])
    }
    ci.PC <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)

    par(mfrow=c(1,1))
    plot(x=list(),y=list(),xlim=c(100,500),ylim=c(-0.2,1.2),ylab="Value",xlab="Day of Year",main=siteName,cex.axis=2,cex.lab=2,cex.main=2)
    lines(xseq,ci.GOES[2,],col="black",lwd=2)
    lines(xseq,ci.GOES[1,],col="black", lty = 2,lwd=2)
    lines(xseq,ci.GOES[3,],col="black", lty = 2,lwd=2)
    abline(v=GOES.trans,col="black")
    #points(data.GOES$x,rescale(var.GOES.c,var.GOES.d,data.GOES$y),col="Black",pch=20)

    lines(xseq,ci.PC[2,],col="cyan",lwd=2)
    lines(xseq,ci.PC[1,],col="cyan", lty = 2,lwd=2)
    lines(xseq,ci.PC[3,],col="cyan", lty = 2,lwd=2)
    abline(v=PC.trans,col="cyan")
    #points(data.PC$x,rescale(var.PC.c,var.PC.d,data.PC$y),col="cyan",pch=20)

    lines(xseq,ci.MODIS[2,],col="red",lwd=2)
    lines(xseq,ci.MODIS[1,],col="red", lty = 2,lwd=2)
    lines(xseq,ci.MODIS[3,],col="red", lty = 2,lwd=2)
    abline(v=MODIS.trans,col="red")
    #points(data.MODIS$x,rescale(var.MODIS.c,var.MODIS.d,data.MODIS$y),col="red",pch=20)

  }
  dev.off()
}
