##' Create the credible interval envelope for plotting
##' 
##' @param x time range
##' @param ylo the bottom credible interval values
##' @param yhi the top credible interval values
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}

library("rjags")
library("PhenologyBayesModeling")
library("grDevices")
startDay <- 182
endDay <- 546
startDate <- as.Date(startDay,origin="2016-12-31")
endDate <- as.Date(endDay,origin="2016-12-31")
xseq <- seq(startDay,endDay,1)

pdf(file="GOES_ALL_NEWDiurnalFits_russellSage.pdf",width=6,height=10)
par(mfrow=c(5,1),mai=c(0.35,1,0.15,0.2))#,mai=c(0.4,0.4,0.2,0.2))

##Plot All diurnal fits
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
s <- 9
#for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData[s,1])
  c.vals <- numeric()
  days <- numeric()
  counts <- numeric()
  inFileName <- paste(siteName,"_diurnal5FitData.RData",sep="")
  if(file.exists(inFileName)){
    load(inFileName)
    dat <- matrix(nrow=5,ncol=length(data$x))
    dat[1,] <- data$x
    dat[2,] <- data$y
    dat[3,] <- data$obs.prec
    #dat[4,] <- data$size
    sd <- 1/sqrt(data$obs.prec)
    dat[4,] <- data$y - 1.96 * sd
    dat[5,] <- data$y + 1.96 * sd
    
    dat2 <- dat
    dat.dates <- as.Date(dat2[1,],origin="2016-12-31")
    
    plot(dat.dates,dat2[2,],pch=20,col="black",ylab="Diurnal Fit",cex.axis=2,cex.lab=2,xlab="",ylim=c(0,1.2),main=siteName)
    polygon(x=c(as.Date("2017-11-30"),as.Date("2017-12-11"),as.Date("2017-12-11"),as.Date("2017-11-30")),y=c(-0.2,-0.2,1.4,1.4),col=adjustcolor("#fff7bc",0.5),border=NA)
    
    #plot(dat.dates,dat2[2,],pch=20,col="white",ylab="NDVI",cex.axis=2,cex.lab=2,xlab="",ylim=c(0,1),col.axis="white",col.lab="white",cex.lab=3,cex.axis=3)
    #polygon(x=c(as.Date("2017-11-30"),as.Date("2017-12-11"),as.Date("2017-12-11"),as.Date("2017-11-30")),y=c(-0.2,-0.2,1.4,1.4),col=adjustcolor("#fff7bc",0.5),border=NA)
    
    
    OverallinFileName <- paste(siteName,"_overall5_varBurn.RData",sep="")
    print(OverallinFileName)
    if(file.exists(OverallinFileName)){
      load(OverallinFileName)
      var.mat <- as.matrix(var.burn)
      colnames(var.mat)
      GOES.TranF<-var.mat[,1]
      GOES.TranS<-var.mat[,2]
      GOES.bF <- var.mat[,3]
      GOES.bS <- var.mat[,4]
      GOES.c<-var.mat[,5]
      GOES.d<-var.mat[,6]
      GOES.k <- var.mat[,7]
      GOES.prec <- var.mat[,8]
      
      ci.Overall <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = FALSE)
      ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="#d95f02")
    
      }
    for(i in 1:ncol(dat2)){
      #print(dat2[3,i])
      ys <- seq(dat2[4,i],dat2[5,i],0.0001)
      #ys <- seq((dat2[2,i]-0.5*dat2[3,i]),(dat2[2,i]+0.5*dat2[3,i]),0.0001)
      #ys <- seq((dat2[2,i]-1.96*dat2[3,i]),(dat2[2,i]+1.96*dat2[3,i]),0.0001)
      #lines(rep(dat.dates[i],length(ys)),ys,col=adjustcolor("gray",0.65))
      lines(rep(dat.dates[i],length(ys)),ys,col=adjustcolor("gray",0.2))
    }
    #points(dat.dates,dat2[2,],pch=20,col="black")
    points(dat.dates,dat2[2,],pch=20,col="white")
  }
#}
dev.off()
