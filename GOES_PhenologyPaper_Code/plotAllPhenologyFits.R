#' Plot all the phenology fits
#'
#' @return
#' @export
#' @import grDevices
#' @import PhenologyBayesModeling
#' @import ecoforecastR
#'
#' @examples
plotAllPhenologyFits <- function(year,startDate,endDate,siteData){
  GOES.col <- "#1b9e77" ##Teal
  MODIS.N.col <- "#7570b3" ##purple
  MODIS.E.col <- "#d95f02"##orange
  PC.col <- "#e7298a"##pink
  
  xseq <- seq(1,365)
  pdf(paste("GOES_Phenology_Paper_CompareSourcesUpdated",year,".pdf",sep=""),width=20,height=11)  

  par(mfrow=c(5,5),mai=c(0.30,0.6,0.01,0.2))#,mai=c(0.4,0.4,0.2,0.2))

  for(s in 1:nrow(siteData)){
    siteName <- as.character(siteData$siteName[s])
    print(siteName)
    if(!(siteName=="russellsage"&&yr==2019)){
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
      inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_varBurn.RData",sep="") 
      load(inputFileName)
      
      var.mat<-as.matrix(PC.md.out)
      ci.PC <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale=FALSE,seasonOrder="SF")
      ci.PC2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale=TRUE,seasonOrder="SF")
      PC.x <- data.PC$x#[data.PC$x>181]
      PC.y <- data.PC$y#[data.PC$x>181]
      PC.dates <- as.Date(PC.x,origin=paste(as.character(yr-1),"-12-31",sep=""))
      
      plot(PC.dates,PC.y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main=siteName,bty = 'n',ylim=c(0.3,0.6))
      ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.PC[1,],ci.PC[3,],col=PC.col)#"#1f78b4")
      points(PC.dates,PC.y,pch=20)
      ##Need CIs for points 
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
      points(MODIS.dates,data.MODIS.N$y,pch=20)
      ##Need CIs on points 
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
      points(MODIS.dates,data.MODIS.E$y,pch=20)
      print("Done with MODIS EVI")
      
      plot(filteredData$dates,rep(NA,length(filteredData$dates)),ylim=c(0,1),pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,
           xlab="",main="",bty = 'n')
      
      ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.N2[1,],ci.MODIS.N2[3,],col=MODIS.N.col)
      ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.E2[1,],ci.MODIS.E2[3,],col=adjustcolor(MODIS.E.col,0.6))
      ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall2[1,],ci.Overall2[3,],col=adjustcolor(GOES.col,0.6))
      ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.PC2[1,],ci.PC2[3,],col=adjustcolor(PC.col,0.6))#"#1f78b4")
    }
  }
  dev.off()
}

# PC.TranF<-var.mat[,1]
# PC.TranS<-var.mat[,2]
# PC.bF <- var.mat[,3]
# PC.bS <- var.mat[,4]
# PC.c<-var.mat[,5]
# PC.d<-var.mat[,6]
# PC.k <- var.mat[,7]
# PC.prec <- var.mat[,8]

# colnames(var.mat)
# GOES.TranF<-var.mat[,1]
# GOES.TranS<-var.mat[,2]
# GOES.bF <- var.mat[,3]
# GOES.bS <- var.mat[,4]
# GOES.c<-var.mat[,5]
# GOES.d<-var.mat[,6]
# GOES.k <- var.mat[,7]
# GOES.prec <- var.mat[,8]

# MODIS.N.TranF<-var.mat[,1]
# MODIS.N.TranS<-var.mat[,2]
# MODIS.N.bF <- var.mat[,3]
# MODIS.N.bS <- var.mat[,4]
# MODIS.N.c<-var.mat[,5]
# MODIS.N.d<-var.mat[,6]
# MODIS.N.k <- var.mat[,7]
# MODIS.N.prec <- var.mat[,8]

# MODIS.E.TranF<-var.mat[,1]
# MODIS.E.TranS<-var.mat[,2]
# MODIS.E.bF <- var.mat[,3]
# MODIS.E.bS <- var.mat[,4]
# MODIS.E.c<-var.mat[,5]
# MODIS.E.d<-var.mat[,6]
# MODIS.E.k <- var.mat[,7]
# MODIS.E.prec <- var.mat[,8]
