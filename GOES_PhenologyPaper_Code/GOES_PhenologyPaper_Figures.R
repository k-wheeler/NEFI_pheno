
library("grDevices")

GOES.col <- "#1b9e77" ##Teal
MODIS.N.col <- "#7570b3" ##purple
MODIS.E.col <- "#d95f02"##orange
PC.col <- "#e7298a"##pink

#' Figure 1 Code
#' 
#' @param siteData
#'
#' @return
#' @export
#'
#' @examples
createFigure1 <- function(siteData){
  NLCD <- raster('LandCoverData/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img')
  colnames(siteData) <- c("Name","Latitude","Longitude","URL","PFT")
  sites <- siteData[,c(3,2)]
  coordinates(sites) <- c("Longitude","Latitude")
  proj4string(sites) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  sites_transformed <- spTransform(sites, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  plot(NLCD)
  points (sites_transformed, pch=16, col="red", cex=2)
  points (sites_transformed, col="black", cex=2)
}

#' Create First Derivative Values
#'
#' @param TranF 
#' @param bF 
#' @param TranS 
#' @param bS 
#' @param c 
#' @param d 
#' @param k 
#' @param xseq 
#'
#' @return
#' @export
#'
#' @examples
createFirstDerVals <- function(TranF,bF,TranS,bS,c,d,k,xseq){
  xseqSpring <- xseq[xseq<=k]
  num2 <- bS * c * exp(bS*(xseqSpring-TranS))
  den2 <- (exp(bS*(xseqSpring-TranS))+1)^2
  spring <- -1 *num2/den2
  
  xseqFall <- xseq[xseq>k]
  num1 <- bF * c * exp(bF*(xseqFall-TranF))
  den1 <- (exp(bF*(xseqFall-TranF))+1)^2
  fall <- -1 *num1/den1
  return(c(spring,fall))
}

#' Create Second Derivative Values
#'
#' @param TranF 
#' @param bF 
#' @param TranS 
#' @param bS 
#' @param c 
#' @param d 
#' @param k 
#' @param xseq 
#'
#' @return
#' @export
#'
#' @examples
createSecondDerVals <- function(TranF,bF,TranS,bS,c,d,k,xseq){
  xseqSpring <- xseq[xseq<=k]
  num2 <- bS**2*c*(exp(bS*(xseqSpring-TranS))-1)*exp(bS*(xseqSpring-TranS))
  den2 <- (exp(bS*(xseqSpring-TranS))+1)**3
  spring <- -1 *num2/den2
  
  xseqFall <- xseq[xseq>k]
  num1 <- bF**2*c*(exp(bF*(xseqFall-TranF))-1)*exp(bF*(xseqFall-TranF))
  den1 <- (exp(bF*(xseqFall-TranF))+1)**3
  fall <- -1 *num1/den1
  return(c(spring,fall))
}

#' Create Third Derivative Values
#'
#' @param TranF 
#' @param bF 
#' @param TranS 
#' @param bS 
#' @param c 
#' @param d 
#' @param k 
#' @param xseq 
#'
#' @return
#' @export
#'
#' @examples
createThirdDerVals <- function(TranF,bF,TranS,bS,c,d,k,xseq){
  xseqSpring <- xseq[xseq<=k]
  num2 <- bS**3*c*exp(bS*(xseqSpring-TranS))*(exp(2*bS*(xseqSpring-TranS))-4*exp(bS*(xseqSpring-TranS))+1)
  den2 <- (exp(bS*(xseqSpring-TranS))+1)**4
  spring <- -1 * num2/den2
  
  xseqFall <- xseq[xseq>k]
  num1 <- bF**3*c*exp(bF*(xseqFall-TranF))*(exp(2*bF*(xseqFall-TranF))-4*exp(bF*(xseqFall-TranF))+1)
  den1 <- (exp(bF*(xseqFall-TranF))+1)**4
  fall <- num1/den1
  
  return(c(fall,spring))
}

#' Calculate transition roots
#'
#' @param c 
#' @param b 
#' @param m 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
calTransRoots <- function(c,b,m,d){
  root1 <- (b*m+log(sqrt(3)+2))/b
  root2 <- (b*m+log(2-sqrt(3)))/b
  return(c(root1,root2))
}

#' Create Figure 2
#'
#' @return
#' @export
#' @import PhenologyBayesModeling
createFigure2 <- function(){
  xseqDates <- seq(as.Date("2018-01-01"),as.Date("2018-12-31"),"day")
  xseq <- seq(1,365)
  yvals <- deciduousYvals(TranF=300,bF=0.15,TranS=(475-365),bS=-0.2,c=0.2,d=0.2,k=182,xseq=xseq,seasonOrder = "SF")
  yvals1st <- createFirstDerVals(TranF=300,bF=0.15,TranS=(475-365),bS=-0.2,c=0.2,d=0.2,k=182,xseq=xseq)
  yvals2nd <- createSecondDerVals(TranF=300,bF=0.15,TranS=(475-365),bS=-0.2,c=0.2,d=0.2,k=182,xseq=xseq)
  yvals3rd <- createThirdDerVals(TranF=300,bF=0.15,TranS=(475-365),bS=-0.2,c=0.2,d=0.2,k=182,xseq=xseq)
  fallRoots <- round(calTransRoots(c=0.2,b=0.15,m=300,d=0.2),digits=0)
  springRoots <- round(calTransRoots(c=0.2,b=-0.20,m=(475-365),d=0.2),digits=0)
  jpeg("GOES_Phenology_Paper_ExplainTransitionDates_RAW.jpeg",width=5,height=8,units = "in",res=1000)  
  
  par(mfrow=c(4,1),mai=c(0.4,0.1,0.01,0.01))
  plot(xseqDates,yvals,pch=".",bty = 'n')
  #abline(h=0,col="gray",lwd=2)
  lines(xseqDates,yvals,lwd=3,col="gray")
  abline(v=xseqDates[xseq==fallRoots[1]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==fallRoots[2]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==springRoots[1]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==springRoots[2]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==300],lwd=2,lty=1)#,col="purple",lwd=2)
  abline(v=xseqDates[xseq==(475-365)],lwd=2,lty=1)#,col="purple",lwd=2)
  #legend("bottomleft",lty=c(1,1,2),col=c("gray","black","black"),c("Curve","Transition Dates","Value of 0"),cex=2)
  
  plot(xseqDates,yvals1st,pch=".",bty = 'n')
  lines(xseqDates,yvals1st,lwd=3,col="gray")
  abline(h=0,lwd=1,lty=2)
  abline(v=xseqDates[xseq==fallRoots[1]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==fallRoots[2]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==springRoots[1]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==springRoots[2]],lwd=2,lty=1)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==300],lwd=2,lty=1)#col="purple",lwd=2)
  abline(v=xseqDates[xseq==(475-365)],lwd=2,lty=1)#,col="purple",lwd=2)
  
  plot(xseqDates,yvals2nd,pch=".",bty = 'n')
  lines(xseqDates,yvals2nd,lwd=3,col="gray")
  abline(h=0,lwd=1,lty=2)
  abline(v=xseqDates[xseq==fallRoots[1]],lwd=2)
  abline(v=xseqDates[xseq==fallRoots[2]],lwd=2)
  abline(v=xseqDates[xseq==springRoots[1]],lwd=2)
  abline(v=xseqDates[xseq==springRoots[2]],lwd=2)
  abline(v=xseqDates[xseq==300],lwd=2)
  abline(v=xseqDates[xseq==(475-365)],lwd=2)
  
  plot(xseqDates,yvals3rd,pch=".",bty = 'n')
  lines(xseqDates,yvals3rd,lwd=3,col="gray")
  abline(h=0,lwd=1,lty=2)
  abline(v=xseqDates[xseq==fallRoots[1]],lwd=2)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==fallRoots[2]],lwd=2)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==springRoots[1]],lwd=2)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==springRoots[2]],lwd=2)#,col="cyan",lwd=3)
  abline(v=xseqDates[xseq==300],lwd=2)#,col="purple",lwd=2)
  abline(v=xseqDates[xseq==(475-365)],lwd=2)#,col="purple",lwd=2)
  
  dev.off()
  
}

#' Create Figure 3: Example Fits
#'
#' @return
#' @export
#' @import rjags
#' @import PhenologyBayesModeling
#' @import grDevices
#' @import ecoforecastR
#'
#' @examples
createFigure3 <- function(siteData){
  GOES.col <- "#1b9e77" ##Teal
  MODIS.N.col <- "#7570b3" ##purple
  MODIS.E.col <- "#d95f02"##orange
  PC.col <- "#e7298a"##pink
  
  startDate <- as.Date("2018-01-01")
  endDate <- as.Date("2018-12-31")
  xseq <- seq(1,365)
  ####################

  jpeg("GOES_Phenology_Paper_CompareSourcesSelectUpdated_RAW.jpeg",width=9,height=10,units = "in",res=1000)  
  par(mfcol=c(5,2),mai=c(0.30,0.6,0.2,0.2))#,mai=c(0.4,0.4,0.2,0.2))
  
  iseq <- c(12,1)
  years <- c(2018,2019)
  
  for(ival in 1:length(iseq)){
    yr <- years[ival]
    s <- iseq[ival]
    siteName <- as.character(siteData$siteName[s])
    print(siteName)
    URL <- as.character(siteData$URL[s])
    
    #GOES
    OverallinFileName <- paste(siteName,"_",startDate,"_",endDate,"_GOES_varBurn.RData",sep="")
    print(OverallinFileName)
    load(OverallinFileName)
    var.mat <- as.matrix(var.burn)
    
    ci.Overall <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = FALSE,seasonOrder="SF")
    ci.Overall2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    load(file=paste(savePath,siteName,"_",yr,"_diurnalFitDataFiltered.RData",sep=""))
    
    filteredData$dates <- as.Date(filteredData$x,origin=paste(as.character(yr-1),"-12-31",sep=""))
    filteredData$Q1 <- filteredData$y - 1.96 * sqrt(1/filteredData$obs.prec)
    filteredData$Q3 <- filteredData$y + 1.96 * sqrt(1/filteredData$obs.prec)
    
    plot(filteredData$dates,filteredData$y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0.1,0.9))
    for(j in 1:length(filteredData$x)){
      lines(rep(filteredData$dates[j],2),c(filteredData$Q1[j],filteredData$Q3[j]),col=adjustcolor("black",0.3),lwd=1)
    }
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
    points(filteredData$dates,filteredData$y,pch=20,col="black")
    
    print("Done with GOES")
    
    ##Plot PhenoCam
    data.PC = PC_data(siteName=siteName,URL=URL,startDate=startDate,endDate=endDate,seasonOrder="SF")
    data.PC$sd <- sqrt(1/data.PC$obs.prec)
    
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(PC.md.out)
    ci.PC <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale=FALSE,seasonOrder="SF")
    ci.PC2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale=TRUE,seasonOrder="SF")
    PC.x <- data.PC$x
    PC.y <- data.PC$y
    PC.sd <- data.PC$sd
    PC.dates <- as.Date(PC.x,origin=paste(as.character(yr-1),"-12-31",sep=""))
    
    plot(PC.dates,PC.y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0.25,0.5))
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.PC[1,],ci.PC[3,],col=PC.col)#"#1f78b4")
    for(i in 1:length(data.PC$sd)){
      lines(c(PC.dates[i],PC.dates[i]),c((PC.y[i]-(1.96*PC.sd[i]))
                                         ,(PC.y[i]+(1.96*PC.sd[i]))),col=adjustcolor("black",0.3))
    }
    points(PC.dates,PC.y,pch=20)
    
    print("Done with PC")
    
    ##MODIS NDVI
    data.MODIS.N = MODIS_data(siteName=siteName,lat=lat,long=long,startDate=startDate,endDate=endDate,
                              metric="NDVI",seasonOrder="SF",lastYear=yr)
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_NDVI_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.N.md.out)
    ci.MODIS.N <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = FALSE,seasonOrder="SF")
    ci.MODIS.N2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    MODIS.dates <- as.Date(data.MODIS.N$x,origin=paste(as.character(yr-1),"-12-31",sep=""))
    plot(MODIS.dates,data.MODIS.N$y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0.3,1))
    
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.N[1,],ci.MODIS.N[3,],col=MODIS.N.col)
    for(i in 1:length(data.MODIS.N$x)){
      lines(c(MODIS.dates[i],MODIS.dates[i]),c((data.MODIS.N$y[i] - 0.01)
                                               ,(data.MODIS.N$y[i] + 0.01)),col=adjustcolor("black",0.6))
    }
    points(MODIS.dates,data.MODIS.N$y,pch=20)
    print("Done with MODIS NDVI")
    
    ##MODIS EVI
    data.MODIS.E = MODIS_data(siteName=siteName,lat=lat,long=long,startDate=startDate,endDate=endDate,
                              metric="EVI",seasonOrder="SF",lastYear=yr)
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_EVI_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.E.md.out)
    ci.MODIS.E <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = FALSE,seasonOrder="SF")
    ci.MODIS.E2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    MODIS.dates <- as.Date(data.MODIS.E$x,origin=paste(as.character(yr-1),"-12-31",sep=""))
    plot(MODIS.dates,data.MODIS.E$y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0,0.8))
    
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.E[1,],ci.MODIS.E[3,],col=MODIS.E.col)
    for(i in 1:length(data.MODIS.E$x)){
      lines(c(MODIS.dates[i],MODIS.dates[i]),c((data.MODIS.E$y[i] - 0.02)
                                               ,(data.MODIS.E$y[i] + 0.02)),col=adjustcolor("black",0.6))
    }
    points(MODIS.dates,data.MODIS.E$y,pch=20)
    print("Done with MODIS EVI")
    
    plot(filteredData$dates,rep(NA,length(filteredData$dates)),ylim=c(0,1),pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,
         xlab="",main="",bty = 'n')
    
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.N2[1,],ci.MODIS.N2[3,],col=MODIS.N.col)
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.E2[1,],ci.MODIS.E2[3,],col=adjustcolor(MODIS.E.col,0.6))
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall2[1,],ci.Overall2[3,],col=adjustcolor(GOES.col,0.6))
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.PC2[1,],ci.PC2[3,],col=adjustcolor(PC.col,0.6))#"#1f78b4")
  }
  dev.off()
  
}

plotCIs <- function(CIs,ylim,label){
  #colnames(CIs) <- c("CI","model")
  #plot(as.character(CIs.TranSs[,2]),as.numeric(CIs.TranSs[,1]))
  
  # CI <- CIs[CIs$model=="GOES",1]
  # plot(density(as.numeric(CI),from=0),col=GOES.col,lwd=2,ylim=c(0,0.045))
  # CI <- CIs[CIs$model=="PC",1]
  # lines(density(as.numeric(CI),from=0),col=PC.col,lwd=2)
  # CI <- CIs[CIs$model=="MODISN",1]
  # lines(density(as.numeric(CI),from=0),col=MODIS.N.col,lwd=2)
  # CI <- CIs[CIs$model=="MODISE",1]
  # lines(density(as.numeric(CI),from=0),col=MODIS.E.col,lwd=2)
  
  #plot(CIs$model,as.numeric(CIs$CI))
  boxplot(as.numeric(CI)~model, data=CIs,pch=20,col=PC.col,ylim=ylim,frame.plot=FALSE,xaxt="n",yaxt="n",cex=1.5)
  boxplot(as.numeric(CI)~model, data=CIs[CIs$model=="GOES",],pch=20,add=TRUE,col=GOES.col,frame.plot=FALSE,xaxt="n",yaxt="n",cex=1.5)
  boxplot(as.numeric(CI)~model, data=CIs[CIs$model=="MODISN",],pch=20,add=TRUE,col=MODIS.N.col,frame.plot=FALSE,xaxt="n",yaxt="n",cex=1.5)
  boxplot(as.numeric(CI)~model, data=CIs[CIs$model=="MODISE",],pch=20,add=TRUE,col=MODIS.E.col,frame.plot=FALSE,xaxt="n",yaxt="n",cex=1.5)
  
  axisVals <- seq(ylim[1],ylim[2],20)
  print(axisVals)
  axis(2,at=axisVals,labels = FALSE,cex.axis=2)
}

#' Create Figure 4
#'
#' @return
#' @export
#'
#' @examples
createFigure4 <- function(){
  ##Calculate the transition date CIs 
  ##Spring
  GOES.col <- "#1b9e77" ##Teal
  MODIS.N.col <- "#7570b3" ##purple
  MODIS.E.col <- "#d95f02"##orange
  PC.col <- "#e7298a"##pink
  
  GOES.CIs.TranSs1 <- cbind(data.frame(apply(GOES.TranSs1,1,quantile,0.975)-apply(GOES.TranSs1,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranSs1) <- c("CI","model")
  GOES.CIs.TranSs <- cbind(data.frame(apply(GOES.TranSs,1,quantile,0.975)-apply(GOES.TranSs,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranSs) <- c("CI","model")
  GOES.CIs.TranSs2 <- cbind(data.frame(apply(GOES.TranSs2,1,quantile,0.975)-apply(GOES.TranSs2,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranSs2) <- c("CI","model")
  
  PC.CIs.TranSs1 <- cbind(data.frame(apply(PC.TranSs1,1,quantile,0.975)-apply(PC.TranSs1,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranSs1) <- c("CI","model")
  PC.CIs.TranSs <- cbind(data.frame(apply(PC.TranSs,1,quantile,0.975)-apply(PC.TranSs,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranSs) <- c("CI","model")
  PC.CIs.TranSs2 <- cbind(data.frame(apply(PC.TranSs2,1,quantile,0.975)-apply(PC.TranSs2,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranSs2) <- c("CI","model")
  
  MODIS.E.CIs.TranSs1 <- cbind(data.frame(apply(MODIS.E.TranSs1,1,quantile,0.975)-apply(MODIS.E.TranSs1,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranSs1) <- c("CI","model")
  MODIS.E.CIs.TranSs <- cbind(data.frame(apply(MODIS.E.TranSs,1,quantile,0.975)-apply(MODIS.E.TranSs,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranSs) <- c("CI","model")
  MODIS.E.CIs.TranSs2 <- cbind(data.frame(apply(MODIS.E.TranSs2,1,quantile,0.975)-apply(MODIS.E.TranSs2,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranSs2) <- c("CI","model")
  
  MODIS.N.CIs.TranSs1 <- cbind(data.frame(apply(MODIS.N.TranSs1,1,quantile,0.975)-apply(MODIS.N.TranSs1,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranSs1) <- c("CI","model")
  MODIS.N.CIs.TranSs <- cbind(data.frame(apply(MODIS.N.TranSs,1,quantile,0.975)-apply(MODIS.N.TranSs,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranSs) <- c("CI","model")
  MODIS.N.CIs.TranSs2 <- cbind(data.frame(apply(MODIS.N.TranSs2,1,quantile,0.975)-apply(MODIS.N.TranSs2,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranSs2) <- c("CI","model")
  
  GOES.CIs.TranFs1 <- cbind(data.frame(apply(GOES.TranFs1,1,quantile,0.975)-apply(GOES.TranFs1,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranFs1) <- c("CI","model")
  GOES.CIs.TranFs <- cbind(data.frame(apply(GOES.TranFs,1,quantile,0.975)-apply(GOES.TranFs,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranFs) <- c("CI","model")
  GOES.CIs.TranFs2 <- cbind(data.frame(apply(GOES.TranFs2,1,quantile,0.975)-apply(GOES.TranFs2,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranFs2) <- c("CI","model")
  
  PC.CIs.TranFs1 <- cbind(data.frame(apply(PC.TranFs1,1,quantile,0.975)-apply(PC.TranFs1,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranFs1) <- c("CI","model")
  PC.CIs.TranFs <- cbind(data.frame(apply(PC.TranFs,1,quantile,0.975)-apply(PC.TranFs,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranFs) <- c("CI","model")
  PC.CIs.TranFs2 <- cbind(data.frame(apply(PC.TranFs2,1,quantile,0.975)-apply(PC.TranFs2,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranFs2) <- c("CI","model")
  
  MODIS.N.CIs.TranFs1 <- cbind(data.frame(apply(MODIS.N.TranFs1,1,quantile,0.975)-apply(MODIS.N.TranFs1,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranFs1) <- c("CI","model")
  MODIS.N.CIs.TranFs <- cbind(data.frame(apply(MODIS.N.TranFs,1,quantile,0.975)-apply(MODIS.N.TranFs,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranFs) <- c("CI","model")
  MODIS.N.CIs.TranFs2 <- cbind(data.frame(apply(MODIS.N.TranFs2,1,quantile,0.975)-apply(MODIS.N.TranFs2,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranFs2) <- c("CI","model")
  
  MODIS.E.CIs.TranFs1 <- cbind(data.frame(apply(MODIS.E.TranFs1,1,quantile,0.975)-apply(MODIS.E.TranFs1,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranFs1) <- c("CI","model")
  MODIS.E.CIs.TranFs <- cbind(data.frame(apply(MODIS.E.TranFs,1,quantile,0.975)-apply(MODIS.E.TranFs,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranFs) <- c("CI","model")
  MODIS.E.CIs.TranFs2 <- cbind(data.frame(apply(MODIS.E.TranFs2,1,quantile,0.975)-apply(MODIS.E.TranFs2,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranFs2) <- c("CI","model")
  
  CIs <- data.frame(rbind(PC.CIs.TranSs,GOES.CIs.TranSs,MODIS.N.CIs.TranSs,MODIS.E.CIs.TranSs))
  
  
  jpeg("GOES_Phenology_Paper_CIwidths_RAW1Updated.jpeg",width=3,height=2,units = "in",res=1000)  
  
  par(mfrow=c(1,3),mai=c(0.01,0.1,0.1,0.01))
  plotCIs(CIs = data.frame(rbind(PC.CIs.TranSs1,GOES.CIs.TranSs1,MODIS.N.CIs.TranSs1,MODIS.E.CIs.TranSs1)),ylim=c(0,60),label=TRUE)
  plotCIs(CIs = data.frame(rbind(PC.CIs.TranSs,GOES.CIs.TranSs,MODIS.N.CIs.TranSs,MODIS.E.CIs.TranSs)),ylim=c(0,60),label=FALSE)
  plotCIs(CIs = data.frame(rbind(PC.CIs.TranSs2,GOES.CIs.TranSs2,MODIS.N.CIs.TranSs2,MODIS.E.CIs.TranSs2)),ylim=c(0,60),label=FALSE)
  dev.off()
  
  jpeg("GOES_Phenology_Paper_CIwidths_RAW2Updated.jpeg",width=3,height=5,units = "in",res=1000)  
  
  par(mfrow=c(1,3),mai=c(0.01,0.1,0.1,0.01))
  
  plotCIs(CIs = data.frame(rbind(PC.CIs.TranFs1,GOES.CIs.TranFs1,MODIS.N.CIs.TranFs1,MODIS.E.CIs.TranFs1)),ylim=c(0,150),label=TRUE)
  plotCIs(CIs = data.frame(rbind(PC.CIs.TranFs,GOES.CIs.TranFs,MODIS.N.CIs.TranFs,MODIS.E.CIs.TranFs)),ylim=c(0,150),label=FALSE)
  plotCIs(CIs = data.frame(rbind(PC.CIs.TranFs2,GOES.CIs.TranFs2,MODIS.N.CIs.TranFs2,MODIS.E.CIs.TranFs2)),ylim=c(0,150),label=FALSE)
  
  dev.off()
}

#' Create Figure 5
#'
#' @return
#' @export
#'
createFigure5 <- function(){
  ##Figure 5: Transition date figure 
  ######Note: transition date matrices created in the calculateCompareMetrics.R file
  
  GOES.col <- adjustcolor("#1b9e77",0.6) ##Teal
  MODIS.N.col <- adjustcolor("#7570b3",0.7) ##purple
  MODIS.E.col <- adjustcolor("#d95f02",0.6)##orange
  PC.col <- adjustcolor("#e7298a",0.6)##pink
  
  jpeg("GOES_Phenology_Paper_Transition_Figure_RAWUpdated.jpeg",width=9,height=6,units = "in",res=1000)  
  par(mfrow=c(2,3),pty="s",mai = c(0.3, 0.5, 0.12, 0.12))
  ##Spring Start
  plot(c(),c(),ylim=c(60,140),xlim=c(60,140),ylab="",cex.axis=2.5,xlab="",frame.plot = FALSE)
  abline(a=0,b=1,col="black")

  xMatrixALL <- PC.TranSs1
  yMatrixALL <- MODIS.N.TranSs1
  curCol <- MODIS.N.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  xMatrixALL <- PC.TranSs1
  yMatrixALL <- MODIS.E.TranSs1
  curCol <- MODIS.E.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  xMatrixALL <- PC.TranSs1
  yMatrixALL <- GOES.TranSs1
  curCol <- GOES.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  ##Spring Middle
  plot(c(),c(),ylim=c(60,160),xlim=c(60,160),ylab="",cex.axis=2.5,xlab="",frame.plot = FALSE)
  abline(a=0,b=1,col="black")
  
  xMatrixALL <- PC.TranSs
  yMatrixALL <- MODIS.N.TranSs
  curCol <- MODIS.N.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  xMatrixALL <- PC.TranSs
  yMatrixALL <- MODIS.E.TranSs
  curCol <- MODIS.E.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  xMatrixALL <- PC.TranSs
  yMatrixALL <- GOES.TranSs
  curCol <- GOES.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  ##Spring End
  plot(c(),c(),ylim=c(80,160),xlim=c(80,160),ylab="",cex.axis=2.5,xlab="",frame.plot = FALSE)
  abline(a=0,b=1,col="black")
  
  xMatrixALL <- PC.TranSs2
  yMatrixALL <- MODIS.N.TranSs2
  curCol <- MODIS.N.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  xMatrixALL <- PC.TranSs2
  yMatrixALL <- MODIS.E.TranSs2
  curCol <- MODIS.E.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  xMatrixALL <- PC.TranSs2
  yMatrixALL <- GOES.TranSs2
  curCol <- GOES.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  ############Autumn
  ##Autumn Start
  plot(c(),c(),ylim=c(160,350),xlim=c(160,350),ylab="",cex.axis=2.5,xlab="",frame.plot = FALSE)
  abline(a=0,b=1,col="black")
  
  xMatrixALL <- PC.TranFs1
  yMatrixALL <- MODIS.N.TranFs1
  curCol <- MODIS.N.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  xMatrixALL <- PC.TranFs1
  yMatrixALL <- MODIS.E.TranFs1
  curCol <- MODIS.E.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  xMatrixALL <- PC.TranFs1
  yMatrixALL <- GOES.TranFs1
  curCol <- GOES.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  ##Autumn Middle
  plot(c(),c(),ylim=c(220,360),xlim=c(220,360),ylab="",cex.axis=2.5,xlab="",frame.plot = FALSE)
  abline(a=0,b=1,col="black")
  
  xMatrixALL <- PC.TranFs
  yMatrixALL <- MODIS.N.TranFs
  curCol <- MODIS.N.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  xMatrixALL <- PC.TranFs
  yMatrixALL <- MODIS.E.TranFs
  curCol <- MODIS.E.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  xMatrixALL <- PC.TranFs
  yMatrixALL <- GOES.TranFs
  curCol <- GOES.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  ##Autumn End
  plot(c(),c(),ylim=c(250,370),xlim=c(250,370),ylab="",cex.axis=2.5,xlab="",frame.plot = FALSE)
  abline(a=0,b=1,col="black")
  
  xMatrixALL <- PC.TranFs2
  yMatrixALL <- MODIS.N.TranFs2
  curCol <- MODIS.N.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  
  xMatrixALL <- PC.TranFs2
  yMatrixALL <- MODIS.E.TranFs2
  curCol <- MODIS.E.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  #abline(h=365,col="red")
  xMatrixALL <- PC.TranFs2
  yMatrixALL <- GOES.TranFs2
  curCol <- GOES.col
  xMatrix <- t(apply(xMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  yMatrix <- t(apply(yMatrixALL,FUN=quantile,MARGIN=1,c(0.025,0.5,0.975)))
  
  points(xMatrix[,2],yMatrix[,2],col=curCol,pch=20,cex=2)
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,1],xMatrix[i,3]),c(yMatrix[i,2],yMatrix[i,2]),col=adjustcolor(curCol,0.7))
  }
  for(i in 2:nrow(xMatrix)){
    lines(c(xMatrix[i,2],xMatrix[i,2]),c(yMatrix[i,1],yMatrix[i,3]),col=adjustcolor(curCol,0.7))
  }
  dev.off()
}

#' Add Box Plot (to Figure 6 graph)
#'
#' @param tranName 
#' @param ylim 
#' @param label 
#'
#' @return
#' @export
#'
#' @examples
addBoxPlot <- function(tranName,ylim,label){
  GOES.col <- "#1b9e77" ##Teal
  MODIS.N.col <- "#7570b3" ##purple
  MODIS.E.col <- "#d95f02"##orange
  PC.col <- "#e7298a"##pink
  dat <- read.csv(paste("PhenologyBias_GOES_PC_",tranName,".csv",sep=""),header=TRUE)
  dat <- cbind(dat,rep("GOES_PC",nrow(dat)))
  colnames(dat) <- c("siteName","Q1","Q2","Q3","model")
  subDat <- read.csv(paste("PhenologyBias_MODISN_PC_",tranName,".csv",sep=""),header=TRUE)
  subDat <- cbind(subDat,rep("MODISN_PC",nrow(subDat)))
  colnames(subDat) <- c("siteName","Q1","Q2","Q3","model")

  dat <- rbind(dat,subDat)
  subDat <- read.csv(paste("PhenologyBias_MODISE_PC_",tranName,".csv",sep=""),header=TRUE)
  subDat <- cbind(subDat,rep("MODISE_PC",nrow(subDat)))
  colnames(subDat) <- c("siteName","Q1","Q2","Q3","model")
  dat <- rbind(dat,subDat)

  plot(Q2~model,data=dat,col=GOES.col,pch=20,ylim=ylim,frame.plot=FALSE,xaxt="n",yaxt="n")
  boxplot(Q2~model, data=dat[dat$model=="MODISN_PC",],add=TRUE,col=MODIS.N.col,frame.plot=FALSE,xaxt="n",yaxt="n")
  boxplot(Q2~model, data=dat[dat$model=="MODISE_PC",],add=TRUE,col=MODIS.E.col,frame.plot=FALSE,xaxt="n",yaxt="n")

  axisVals <- seq(ylim[1]+10,ylim[2],20)
  axis(2,at=axisVals,labels = label,cex.axis=2)
  abline(h=0,col="red",lwd=3)
}

#' Create Figure 6
#'
#' @return
#' @export
#'
#' @examples
createFigure6 <- function(){
  #jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_biasPlotsAutumn_RAW.jpeg",width=6,height=2.5,units = "in",res=1000)  
  jpeg("GOES_Phenology_Paper_biasPlots_RAW1Updated.jpeg",width=8,height=2.5,units = "in",res=1000)  
  
  par(mfrow=c(1,3),mai=c(0.01,0.4,0.01,0.1))
  addBoxPlot(tranName="springTran1",ylim=c(-30,20),label=TRUE) #Second addBoxPlot() function
  addBoxPlot(tranName="spring50",ylim=c(-30,20),label=FALSE)
  addBoxPlot(tranName="springTran2",ylim=c(-30,20),label=FALSE)
  
  dev.off()
  
  jpeg("GOES_Phenology_Paper_biasPlots_RAW2Updated.jpeg",width=8,height=9,units = "in",res=1000)  
  
  par(mfrow=c(1,3),mai=c(0.01,0.4,0.01,0.1))
  addBoxPlot(tranName="fallTran1",ylim=c(-90,80),label=TRUE)
  addBoxPlot(tranName="fall50",ylim=c(-90,80),label=FALSE)
  addBoxPlot(tranName="fallTran2",ylim=c(-90,80),label=FALSE)
  
  dev.off()
}
