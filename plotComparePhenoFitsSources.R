library("grDevices")
library("rjags")
library("PhenologyBayesModeling")

##' Create the credible interval envelope for plotting
##' 
##' @param x time range
##' @param ylo the bottom credible interval values
##' @param yhi the top credible interval values
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}



pdf(file="GOES_Pheno_CompareSources4.pdf",width=12,height=8)
par(mfrow=c(2,2),mai=c(0.35,1,0.3,0.2))#,mai=c(0.4,0.4,0.2,0.2))
#sites <- c("HarvardForest","russellSage","shiningrock","luckyHills")
#PFTs <- c("DB","DB","DB","DB")
#sites <- c("HarvardForest","russellSage","shiningrock","greenridge")
sites <- as.character(read.csv("GOES_Paper_sites.csv",header=TRUE)[c(seq(1,6),seq(8,10),seq(15,20)),1])
PFTs <- as.character(read.csv("GOES_Paper_sites.csv",header=TRUE)[c(seq(1,6),seq(8,10),seq(15,20)),5])
for(s in 1:length(sites)){
  PFT=PFTs[s]
  siteName <- sites[s]
  if(siteName=="HarvardForest"){
    siteTitle <- "Harvard Forest"
  }
  else if(siteName=="russellSage"){
    siteTitle <- "Russell Sage"
  }
  else if(siteName=="shiningrock"){
    siteTitle <- "Shining Rock"
  }
  else if(siteName=="greenridge"){
    siteTitle <- "Green Ridge"
  }
  
  if(PFTs[s]=="DB"){
    startDay <- 182
    endDay <- 546
    startDate <- as.Date(startDay,origin="2016-12-31")
    endDate <- as.Date(endDay,origin="2016-12-31")
    xseq <- seq(startDay,endDay,1)
  }
  else if(PFTs[s]=="SH"){
    startDay <- 110
    endDay <- 455
    startDate <- as.Date(startDay,origin="2016-12-31")
    endDate <- as.Date(endDay,origin="2016-12-31")
    xseq <- seq(startDay,endDay,1)
  }
  
  ##GOES Data
  BasicInFileName <- paste(siteName,"_GOES_varBurn.RData",sep="")
  noonDat <- read.csv(paste("GOES_NDVI_",siteName,"_",startDate,"_",endDate,"_noon.csv",sep=""),header=FALSE)
  #noonDat <- read.csv("GOES_NDVI_HarvardForest_2017-07-01_2018-06-30_noon.csv",header=FALSE)
  for(i in 1:length(noonDat)){
    if(noonDat[1,i]<startDay){
      noonDat[1,i] <- noonDat[1,i]+365
    }
  }
  noon.dates <- as.Date(as.numeric(noonDat[1,]),origin="2016-12-31")
  plot(noon.dates,rep(NA,length(noon.dates)),ylim=c(0,1),ylab="Value",pch=20,col="black",main=siteName,cex.axis=2,cex.lab=2,xlab="",cex.main=2)
  OverallinFileName <- paste(siteName,"_overall6_varBurn.RData",sep="")
  
  print(OverallinFileName)
  load(OverallinFileName)
  var.mat <- as.matrix(var.burn)
  print(colnames(var.mat))
  GOES.TranF<-var.mat[,1]
  GOES.TranS<-var.mat[,2]
  GOES.bF <- var.mat[,3]
  GOES.bS <- var.mat[,4]
  GOES.c<-var.mat[,5]
  GOES.d<-var.mat[,6]
  GOES.k <- var.mat[,7]
  GOES.prec <- var.mat[,8]
  
  ci.Overall <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE)
  #ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="purple")
  print("DONE with GOES")
  ##MODIS
  #data.MODIS.N = MODIS_data(siteName=siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="NDVI")
  inputFileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_NDVI_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-as.matrix(MODIS.N.md.out)
  
  MODIS.N.TranF<-var.mat[,1]
  MODIS.N.TranS<-var.mat[,2]
  MODIS.N.bF <- var.mat[,3]
  MODIS.N.bS <- var.mat[,4]
  MODIS.N.c<-var.mat[,5]
  MODIS.N.d<-var.mat[,6]
  MODIS.N.k <- var.mat[,7]
  MODIS.N.prec <- var.mat[,8]
  ci.MODIS.N <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE)
  
  #MODIS.dates <- as.Date(data.MODIS.N$x,origin="2016-12-31")
  #ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="#a6cee3")
  #plot(MODIS.dates,data.MODIS.N$y,pch=20,col="black",ylab="MODIS",cex.axis=2,cex.lab=2,xlab="")
  
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.MODIS.N[1,],ci.MODIS.N[3,],col=adjustcolor("#b2df8a",0.5))
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col=adjustcolor("purple",0.5))
  print("DONE with MODIS NDVI")
  
  ##MODIS EVI
  #data.MODIS.E = MODIS_data(siteName=siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="EVI")
  inputFileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_EVI_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-as.matrix(MODIS.E.md.out)
  
  MODIS.E.TranF<-var.mat[,1]
  MODIS.E.TranS<-var.mat[,2]
  MODIS.E.bF <- var.mat[,3]
  MODIS.E.bS <- var.mat[,4]
  MODIS.E.c<-var.mat[,5]
  MODIS.E.d<-var.mat[,6]
  MODIS.E.k <- var.mat[,7]
  MODIS.E.prec <- var.mat[,8]
  ci.MODIS.E <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE)
  
  #MODIS.dates <- as.Date(data.MODIS.E$x,origin="2016-12-31")
  #ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="#a6cee3")
  #plot(MODIS.dates,data.MODIS.N$y,pch=20,col="black",ylab="MODIS",cex.axis=2,cex.lab=2,xlab="")
  
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.MODIS.E[1,],ci.MODIS.E[3,],col=adjustcolor("orange",0.5))

  print("DONE with MODIS EVI")
  
  
  
  
  ##PC
  #data.PC = PC_data(siteName=siteName,URL="",startDay=startDay,endDay=endDay)
  #inputFileName <- paste(siteName,"_PC_varBurn.RData",sep="")
  inputFileName <- paste(siteName,"_",startDay,"_",endDay,"_PC_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-as.matrix(PC.md.out)
  PC.TranF<-var.mat[,1]
  PC.TranS<-var.mat[,2]
  PC.bF <- var.mat[,3]
  PC.bS <- var.mat[,4]
  PC.c<-var.mat[,5]
  PC.d<-var.mat[,6]
  PC.k <- var.mat[,7]
  PC.prec <- var.mat[,8]
  
  ci.PC <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale=TRUE)
  #PC.x <- data.PC$x[data.PC$x>181]
  #PC.y <- data.PC$y[data.PC$x>181]
  #PC.dates <- as.Date(PC.x,origin="2016-12-31")

  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.PC[1,],ci.PC[3,],col=adjustcolor("#1f78b4",0.5))
  legend("bottomleft",c("GOES Diurnal Fit","MODIS","PhenoCam"),pch=c(20,20,20),col=c("purple","#b2df8a","#1f78b4"))
}
dev.off()
  