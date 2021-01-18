


#################Supplementary Figure on Bayes fits
##Supplementary Figures 2-3
library("grDevices")
library(ecoforecastR)
GOES.col <- "#1b9e77" ##Teal
MODIS.N.col <- "#7570b3" ##purple
MODIS.E.col <- "#d95f02"##orange
PC.col <- "#e7298a"##pink

startDate <- as.Date("2018-01-01")
endDate <- as.Date("2018-12-31")
yr <- 2018
xseq <- seq(1,365)
jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_CompareSourcesALL_Supplementary_RAWUpdated2018.jpeg",width=12,height=14.5,units = "in",res=1000)  

#par(mfrow=c(4,1),mai=c(0.35,1,0.15,0.2))#,mai=c(0.4,0.4,0.2,0.2))
par(mfrow=c(8,2),mai=c(0.30,0.8,0.3,0.2))#,mai=c(0.4,0.4,0.2,0.2))
##Plot Overall GOES
#siteName <- "missouriozarks"
#siteName <- "russellSage"
#siteName <- "shenandoah"
#sites <- as.character(read.csv("GOES_Paper_sites.csv",header=TRUE)[iseq,1])
sites <- as.character(read.csv("GOES_Paper_Sites_FINAL.csv",header=TRUE)$siteName)
siteNames <- c("Harvard Forest","Hubbard Brook","UMBS","Coweeta","Bartlett",
               "Missouri Ozarks", "Morgan Monroe", "Russell Sage", "Willow Creek",
               "Bull Shoals", "Duke", "Green Ridge", "Shenandoah", "Marcell", "Shining Rock")
#q=1
for(i in 1:length(sites)){
  siteName <- sites[i]
  
  #OverallinFileName <- paste(siteName,"_overall6_varBurn.RData",sep="")
  OverallinFileName <- paste(siteName,"_",yr,"_GOES_NDVI_overallFilter_varBurn.RData",sep="")
  print(OverallinFileName)
  load(OverallinFileName)
  var.mat <- as.matrix(var.burn)
  
  ci.Overall <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = FALSE,seasonOrder="SF")
  #ci.Overall2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
  
  load(file=paste(siteName,"_",yr,"_diurnalFitDataModelFilterFilter.RData",sep=""))
  
  filteredData$dates <- as.Date(filteredData$x,origin=paste(as.character(yr-1),"-12-31",sep=""))
  
  filteredData$Q1 <- filteredData$y - 1.96 * sqrt(1/filteredData$obs.prec)
  filteredData$Q3 <- filteredData$y + 1.96 * sqrt(1/filteredData$obs.prec)
  
  
  # plot(filteredData$dates,filteredData$y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0.1,0.9))
  # for(j in 1:length(filteredData$x)){
  #   lines(rep(filteredData$dates[j],2),c(filteredData$Q1[j],filteredData$Q3[j]),col=adjustcolor("black",0.3),lwd=1)
  # }
  # ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
  # points(filteredData$dates,filteredData$y,pch=20,col="black")
  
  if(i %in% c(1,3,5,7,9,11,13,15)){
    ylabV=TRUE
  }else{
    ylabV <- FALSE
  }
  
  if(i %in% c(14,15)){
    xlabV=c("Jan","Mar","May","Jul","Sep","Nov","Jan")
  }else{
    xlabV <- FALSE
  }
  
  plot(filteredData$dates,filteredData$y,pch=20,col="black",ylab=siteNames[i],cex.axis=2,cex.lab=2,
       xlab="",main="",bty = 'n',ylim=c(0.1,0.9),yaxt="n",xaxt="n")
  axis(2,at=c(0.2,0.4,0.6,0.8),labels = ylabV,cex.axis=2)
  
  
  axis(1,at=c(as.Date("2018-01-01"),as.Date("2018-03-01"),as.Date("2018-05-01"),as.Date("2018-07-01"),
              as.Date("2018-09-01"),as.Date("2018-11-01"),as.Date("2019-01-01")),labels = xlabV ,cex.axis=2)
  
  
  
  ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
  for(j in 1:length(filteredData$x)){
    lines(rep(filteredData$dates[j],2),c(filteredData$Q1[j],filteredData$Q3[j]),col="darkgray")##=adjustcolor("black",0.3),lwd=1)
  }
  points(filteredData$dates,filteredData$y,pch=20,col="black")
  
  # ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
  # for(i in 1:ncol(dat2)){
  #   lines(c(dat.dates[i],dat.dates[i]),c(dat2[4,i],dat2[5,i]),col=adjustcolor("black",0.3))
  # }
  # points(dat.dates,dat2[2,],pch=20,col="black")
  
  #print("Done with GOES")
  #q=q+1
  
}
dev.off()


startDate <- as.Date("2019-01-01")
endDate <- as.Date("2019-12-31")
yr <- 2019
xseq <- seq(1,365)
jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_CompareSourcesALL_Supplementary_RAWUpdated2019.jpeg",width=12,height=14.5,units = "in",res=1000)  

#par(mfrow=c(4,1),mai=c(0.35,1,0.15,0.2))#,mai=c(0.4,0.4,0.2,0.2))
par(mfrow=c(8,2),mai=c(0.30,0.8,0.3,0.2))#,mai=c(0.4,0.4,0.2,0.2))
##Plot Overall GOES
#siteName <- "missouriozarks"
#siteName <- "russellSage"
#siteName <- "shenandoah"
#sites <- as.character(read.csv("GOES_Paper_sites.csv",header=TRUE)[iseq,1])
sites <- as.character(read.csv("GOES_Paper_Sites_FINAL.csv",header=TRUE)$siteName)
siteNames <- c("Harvard Forest","Hubbard Brook","UMBS","Coweeta","Bartlett",
               "Missouri Ozarks", "Morgan Monroe", "Russell Sage", "Willow Creek",
               "Bull Shoals", "Duke", "Green Ridge", "Shenandoah", "Marcell", "Shining Rock")
#q=1
for(i in 1:length(sites)){
  siteName <- sites[i]
  
  #OverallinFileName <- paste(siteName,"_overall6_varBurn.RData",sep="")
  OverallinFileName <- paste(siteName,"_",yr,"_GOES_NDVI_overallFilter_varBurn.RData",sep="")
  print(OverallinFileName)
  load(OverallinFileName)
  var.mat <- as.matrix(var.burn)
  
  ci.Overall <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = FALSE,seasonOrder="SF")
  #ci.Overall2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
  
  load(file=paste(siteName,"_",yr,"_diurnalFitDataModelFilterFilter.RData",sep=""))
  
  filteredData$dates <- as.Date(filteredData$x,origin=paste(as.character(yr-1),"-12-31",sep=""))
  
  filteredData$Q1 <- filteredData$y - 1.96 * sqrt(1/filteredData$obs.prec)
  filteredData$Q3 <- filteredData$y + 1.96 * sqrt(1/filteredData$obs.prec)
  
  
  # plot(filteredData$dates,filteredData$y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0.1,0.9))
  # for(j in 1:length(filteredData$x)){
  #   lines(rep(filteredData$dates[j],2),c(filteredData$Q1[j],filteredData$Q3[j]),col=adjustcolor("black",0.3),lwd=1)
  # }
  # ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
  # points(filteredData$dates,filteredData$y,pch=20,col="black")
  
  if(i %in% c(1,3,5,7,9,11,13,15)){
    ylabV=TRUE
  }else{
    ylabV <- FALSE
  }
  
  if(i %in% c(14,15)){
    xlabV=c("Jan","Mar","May","Jul","Sep","Nov","Jan")
  }else{
    xlabV <- FALSE
  }
  
  plot(filteredData$dates,filteredData$y,pch=20,col="black",ylab=siteNames[i],cex.axis=2,cex.lab=2,
       xlab="",main="",bty = 'n',ylim=c(0.1,0.9),yaxt="n",xaxt="n")
  axis(2,at=c(0.2,0.4,0.6,0.8),labels = ylabV,cex.axis=2)
  
  
  axis(1,at=c(as.Date("2019-01-01"),as.Date("2019-03-01"),as.Date("2019-05-01"),as.Date("2019-07-01"),
              as.Date("2019-09-01"),as.Date("2019-11-01"),as.Date("2020-01-01")),labels = xlabV ,cex.axis=2)
  
  
  
  ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
  for(j in 1:length(filteredData$x)){
    lines(rep(filteredData$dates[j],2),c(filteredData$Q1[j],filteredData$Q3[j]),col="darkgray")##=adjustcolor("black",0.3),lwd=1)
  }
  points(filteredData$dates,filteredData$y,pch=20,col="black")
  
  # ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
  # for(i in 1:ncol(dat2)){
  #   lines(c(dat.dates[i],dat.dates[i]),c(dat2[4,i],dat2[5,i]),col=adjustcolor("black",0.3))
  # }
  # points(dat.dates,dat2[2,],pch=20,col="black")
  
  #print("Done with GOES")
  #q=q+1
  
}
dev.off()

###Supplementary Figures 4-5

library("grDevices")
library(PhenologyBayesModeling)
library(ecoforecastR)
GOES.col <- "#1b9e77" ##Teal
MODIS.N.col <- "#7570b3" ##purple
MODIS.E.col <- "#d95f02"##orange
PC.col <- "#e7298a"##pink

xseq <- seq(1,365)
siteData <- read.csv("GOES_Paper_Sites_FINAL.csv",header=TRUE)
siteNames <- c("Harvard Forest","Hubbard Brook","UMBS","Coweeta","Bartlett",
               "Missouri Ozarks", "Morgan Monroe", "Russell Sage", "Willow Creek",
               "Bull Shoals", "Duke", "Green Ridge", "Shenandoah", "Marcell", "Shining Rock")
jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_CompareSourcesCIs2018.jpeg",width=12,height=14.5,units = "in",res=1000)  

#par(mfrow=c(4,1),mai=c(0.35,1,0.15,0.2))#,mai=c(0.4,0.4,0.2,0.2))
par(mfrow=c(8,2),mai=c(0.30,0.6,0.2,0.2))#,mai=c(0.4,0.4,0.2,0.2))
startDate <- as.Date("2018-01-01")
endDate <- as.Date("2018-12-31")
yr <- 2018



for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[s])
  print(siteName)
  if(!(siteName=="russellsage"&&yr==2019)){
    #GOES
    OverallinFileName <- paste(siteName,"_",yr,"_GOES_NDVI_overallFilter_varBurn.RData",sep="")
    print(OverallinFileName)
    load(OverallinFileName)
    
    var.mat <- as.matrix(var.burn)
    
    ci.Overall2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    ##Plot PhenoCam
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_EIV_varBurn_SF.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(PC.md.out)
    
    ci.PC2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale=TRUE,seasonOrder="SF")
    
    
    ##MODIS NDVI
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_NDVI_EIV_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.N.md.out)
    if(siteName=="russellsage"){
      var.mat <- data.frame(var.mat)
      var.mat <- var.mat[var.mat$bF>0,]
    }
    ci.MODIS.N2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    ##MODIS EVI
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_EVI_EIV_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.E.md.out)
    ci.MODIS.E2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    if(s %in% c(1,3,5,7,9,11,13,15)){
      ylabV=TRUE
    }else{
      ylabV <- FALSE
    }
    
    if(s %in% c(14,15)){
      xlabV=c("Jan","Mar","May","Jul","Sep","Nov","Jan")
    }else{
      xlabV <- FALSE
    }
    
    plot(seq(startDate,endDate,"day"),rep(NA,length(seq(startDate,endDate,"day"))),pch=20,col="black",ylab=siteNames[s],cex.axis=2,cex.lab=2,
         xlab="",main="",bty = 'n',ylim=c(0,1),yaxt="n",xaxt="n")
    axis(2,at=c(0.2,0.4,0.6,0.8,1),labels = ylabV,cex.axis=2)
    
    
    axis(1,at=c(as.Date("2018-01-01"),as.Date("2018-03-01"),as.Date("2018-05-01"),as.Date("2018-07-01"),
                as.Date("2018-09-01"),as.Date("2018-11-01"),as.Date("2019-01-01")),labels = xlabV ,cex.axis=2)
    
    
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.N2[1,],ci.MODIS.N2[3,],col=MODIS.N.col)
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.E2[1,],ci.MODIS.E2[3,],col=adjustcolor(MODIS.E.col,0.6))
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall2[1,],ci.Overall2[3,],col=adjustcolor(GOES.col,0.6))
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.PC2[1,],ci.PC2[3,],col=adjustcolor(PC.col,0.8))#"#1f78b4")
  }
}
plot(seq(startDate,endDate,"day"),rep(NA,length(seq(startDate,endDate,"day"))),pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,
     xlab="",main="",bty = 'n',ylim=c(0,1),yaxt="n",xaxt="n")
legend("topleft",lty = c(1,1,1,1),col=c(GOES.col,PC.col,MODIS.N.col,MODIS.E.col),
       c("GOES NDVI", "PhenoCam GCC","MODIS NDVI","MODIS EVI"),lwd=c(5,5,5,5),cex=1.8)
dev.off()

jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_CompareSourcesCIs2019.jpeg",width=12,height=14.5,units = "in",res=1000)  

#par(mfrow=c(4,1),mai=c(0.35,1,0.15,0.2))#,mai=c(0.4,0.4,0.2,0.2))
par(mfrow=c(8,2),mai=c(0.30,0.6,0.2,0.2))#,mai=c(0.4,0.4,0.2,0.2))
startDate <- as.Date("2019-01-01")
endDate <- as.Date("2019-12-31")
yr <- 2019


for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[s])
  print(siteName)
  if(!(siteName=="russellsage"&&yr==2019)){
    #GOES
    OverallinFileName <- paste(siteName,"_",yr,"_GOES_NDVI_overallFilter_varBurn.RData",sep="")
    print(OverallinFileName)
    load(OverallinFileName)
    
    var.mat <- as.matrix(var.burn)
    
    ci.Overall2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    ##Plot PhenoCam
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_EIV_varBurn_SF.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(PC.md.out)
    
    ci.PC2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale=TRUE,seasonOrder="SF")
    
    
    ##MODIS NDVI
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_NDVI_EIV_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.N.md.out)
    ci.MODIS.N2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    ##MODIS EVI
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_EVI_EIV_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.E.md.out)
    ci.MODIS.E2 <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq,doRescale = TRUE,seasonOrder="SF")
    
    if(s %in% c(1,3,5,7,9,11,13,15)){
      ylabV=TRUE
    }else{
      ylabV <- FALSE
    }
    
    if(s %in% c(14,15)){
      xlabV=c("Jan","Mar","May","Jul","Sep","Nov","Jan")
    }else{
      xlabV <- FALSE
    }
    
    plot(seq(startDate,endDate,"day"),rep(NA,length(seq(startDate,endDate,"day"))),pch=20,col="black",ylab=siteNames[s],cex.axis=2,cex.lab=2,
         xlab="",main="",bty = 'n',ylim=c(0,1),yaxt="n",xaxt="n")
    axis(2,at=c(0.2,0.4,0.6,0.8,1),labels = ylabV,cex.axis=2)
    
    
    axis(1,at=c(as.Date("2019-01-01"),as.Date("2019-03-01"),as.Date("2019-05-01"),as.Date("2019-07-01"),
                as.Date("2019-09-01"),as.Date("2019-11-01"),as.Date("2020-01-01")),labels = xlabV ,cex.axis=2)
    
    
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.N2[1,],ci.MODIS.N2[3,],col=MODIS.N.col)
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.MODIS.E2[1,],ci.MODIS.E2[3,],col=adjustcolor(MODIS.E.col,0.6))
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.Overall2[1,],ci.Overall2[3,],col=adjustcolor(GOES.col,0.6))
    ciEnvelope(as.Date(xseq,origin=paste(as.character(yr-1),"-12-31",sep="")),ci.PC2[1,],ci.PC2[3,],col=adjustcolor(PC.col,0.8))#"#1f78b4")
  }
}
plot(seq(startDate,endDate,"day"),rep(NA,length(seq(startDate,endDate,"day"))),pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,
     xlab="",main="",bty = 'n',ylim=c(0,1),yaxt="n",xaxt="n")
plot(seq(startDate,endDate,"day"),rep(NA,length(seq(startDate,endDate,"day"))),pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,
     xlab="",main="",bty = 'n',ylim=c(0,1),yaxt="n",xaxt="n")
legend("topleft",lty = c(1,1,1,1),col=c(GOES.col,PC.col,MODIS.N.col,MODIS.E.col),
       c("GOES NDVI", "PhenoCam GCC","MODIS NDVI","MODIS EVI"),lwd=c(5,5,5,5),cex=1.8)
dev.off()

##Supplementary Figure showing noise filtered out
files <- c("PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnalshiningrock_2019317.csv",
           "PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnalshenandoah_2019061.csv",
           "PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnalmissouriozarks_2019051.csv",
           "PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnalmissouriozarks_2019313.csv",
           "PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnalmissouriozarks_2019314.csv",
           "PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnalmissouriozarks_2019328.csv")

jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_06040Removal.jpeg",width=9.2,height=6.6,units = "in",res=1000)  

par(mfrow=c(2,3),mai=c(0.4,0.4,0.4,0.1))
titles <- c("Shining Rock 13 Nov 2019","Shenandoah 2 Mar 2019","Missouri Ozarks 20 Feb 2019",
            "Missouri Ozarks 9 Nov 2019", "Missouri Ozarks 10 Nov 2019", "Missouri Ozarks 24 Nov 2019")
for(f in 1:length(files)){
  dayData <- read.csv(files[f],header=FALSE)
  plot(as.numeric(dayData[3,]),as.numeric(dayData[2,]),ylim=c(0,1),xlim=c(0,25),pch=20,ylab="NDVI",xlab="Hour",bty="n",main=titles[f])
  abline(h=0.6040,col="red")
  points(as.numeric(dayData[3,]),as.numeric(dayData[2,]),pch=20)
}
dev.off()

dayData2 <- dayData[,!is.na(dayData[2,])]
times <- dayData2[3,round(dayData2[2,],digits=4)==0.6040]
for(t in 1:length(times)){
  hr <- floor(times[t])
  mnt <- (times[t]-floor(times[t]))*60
  paste(as.character(hr),as.character(mnt),sep="")
}
20193170611



##Figure 1: Site Map
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
#library(nneo)
library(googlesheets)
library(ggrepel)
devtools::install_github("ropensci/FedData")
#install.packages("ggrepel")
install.packages("FedData")
library("FedData")
library("raster")

NLCD <- raster('LandCoverData/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img')

# Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
# http://village.anth.wsu.edu
# vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
#                                   proj4string='+proj=utm +datum=NAD83 +zone=12')
# 
# 
# # Get the NLCD (USA ONLY)
# # Returns a raster
# NLCD <- get_nlcd(template=vepPolygon, label='VEPIIN')
# 
# # Plot with raster::plot
# plot(NLCD)


# dta <-  map_data("world",c("usa","canada","mexico"),xlim=c(-100,-65), ylim=c(19,50))

siteData <- data.frame(read.csv("GOES_Paper_Sites.csv",header=TRUE))#,skip=1))
colnames(siteData) <- c("Name","Latitude","Longitude","URL","PFT")
DB.seq <- c(seq(1,6),seq(8,10),seq(15,20))
#SH.seq <- c(7,12,13,14)
siteData <- siteData[DB.seq,]
sites <- siteData[,c(3,2)]
coordinates(sites) <- c("Longitude","Latitude")
proj4string(sites) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
sites_transformed<-spTransform(sites, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#jpeg("GOES_PaperWriting/GOES_Phenology_Paper/GOES_Phenology_Paper_UpdatedMap.jpeg",width=7.6,height=5.96,units="in",res=1000)
plot(NLCD)
points (sites_transformed, pch=16, col="red", cex=2)
points (sites_transformed, col="black", cex=2)
#lbls <- c("Harvard Forest", "Hubbard Brook", "UMBS", "Coweeta", "Bartlett", "Missouri Ozarks","Morgan Monroe","Russell Sage","Willow Creek","Bull Shoals","Duke", "Green Ridge","Shenandoah","Marcell","Shining Rock")


# DB.seq <- c(seq(1,6),seq(8,10),seq(15,20))
# #SH.seq <- c(7,12,13,14)
# siteDataDB <- siteData[DB.seq,]
# siteDataDB$ActualName <- c("Harvard Forest", "Hubbard Brook", "UMBS", "Coweeta", "Bartlett", "Missouri Ozarks","Morgan Monroe","Russell Sage","Willow Creek","Bull Shoals","Duke", "Green Ridge","Shenandoah","Marcell","Shining Rock")
# 
# 
# #dta2 <- dta
# #pdf(file="GOES_Diurnal_SiteLocations.pdf",width=9,height=6.5)
# jpeg("GOES_PaperWriting/GOES_Phenology_Paper/GOES_Phenology_Paper_UpdatedMap.jpeg", width = 6, height = 3.5, units = 'in', res = 1000)
# 
# dta2 <- dta[as.numeric(dta$long)>(-105),]
# #dta2 <- dta2[as.numeric(dta2$long)<(-65),]
# #dta2 <- dta2[as.numeric(dta2$lat)>(20),]
# #dta2 <- dta2[as.numeric(dta2$lat)<(60),]
# gg1 <- ggplot(siteDataDB,main="Selected Site Locations",aes(x=Longitude,y=Latitude))+geom_polygon(data = dta2, aes(x=long, y=lat, group=group),fill=NA,colour="black")# + cord_fixed(1.3)
# gg1 <- gg1 + geom_point(data=siteDataDB, aes(x = Longitude, y = Latitude),color="black",cex=2.5)+theme_bw()
# gg1 <- gg1 + geom_point(data=siteDataDB, aes(x = Longitude, y = Latitude),color="black",cex=2.5,pch=21)+theme_bw()
# #gg1 <- gg1 + geom_point(data=siteData[SH.seq,], aes(x = Longitude, y = Latitude),color="orange",cex=2.5)+theme_bw()
# #gg1 <- gg1 + geom_point(data=siteData[SH.seq,], aes(x = Longitude, y = Latitude),color="black",cex=2.5,pch=21)+theme_bw()
# #gg1 <- gg1 + geom_text(aes(label=siteName),hjust=0, vjust=0)
# gg1 <- gg1 + scale_x_continuous(limits = c(-105.2,-60))
# 
# gg1 <- gg1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# gg1 <- gg1 + labs(x="Longitude",y="Latitude",cex=4)
# gg1 <- gg1 + theme(axis.text=element_text(size=15))
# gg1 <- gg1 + theme(axis.title = element_text(size=20))
# gg1 <- gg1 + theme(panel.background = element_blank())
# 
# gg1 <- gg1 + 
#   geom_label_repel(aes(label = ActualName),
#                    box.padding   = 0.35, 
#                    point.padding = 0.5,
#                    segment.color = 'grey50') +
#   theme_classic()
# xs <- c(-105.2,-104.8,-104.8,-105.2)
# ys <- c(25,25,50,50)
# #blkdt <- list(xs=xs,ys=ys)
# blkdt <- as.data.frame(cbind(xs,ys))
# gg1 <- gg1 + geom_polygon(data=blkdt,aes(x=xs,y=ys),fill="white",colour="white")
# 
# plot(gg1)
# 
# dev.off()





################################################
##Figure 4
##Compare CI widths


########################################################


##Figure 5: p-values


# ##Supplementary Table showing the biases broken down by site
# ##Create sections for different transition dates and have the sites as rows and
#   ##and have different model fits as columns (with quantiles in perenthesis)
# ##What I am picturing is I need to loop over each transition date files 
# 
# combineBiases <- function(pattern){
#   print(pattern)
#   output <- matrix(nrow=16,ncol=0)
#   colNameVals <- "Site Name"
#   
#   sites <- c("Harvard","Hubbard Brook","UMB","Coweeta","Bartlett","Missouri Ozarks","Morgan Monroe",
#              "Russell Sage","Willow Creek","Bull Shoals","Duke","Green Ridge",
#              "Shenandoah", "Marcell", "Shining Rock","Average")
#   output <- cbind(output,sites)
#   
#   biasFiles <- dir(pattern=pattern)
#   for(f in 1:length(biasFiles)){
#     dat <- read.csv(biasFiles[f],header=TRUE)
#     newVal <- paste(round(dat[,3],digits=2)," (", round(dat[,2],digits=2),", ",round(dat[,4],digits=2),")",sep="")
#     output <- cbind(output,newVal)
#     colVal <- paste(strsplit(biasFiles[f],"_")[[1]][2],"_",strsplit(biasFiles[f],"_")[[1]][3],sep="")
#     colNameVals <- c(colNameVals,colVal)
#     colnames(output) <- colNameVals
#   }
#   return(output)
# }
# transDates <- c("springTran1.csv","spring50.csv","springTran2.csv",
#                 "fallTran2.csv","fall50.csv","fallTran1.csv")
# overallOutput <- matrix(nrow=0,ncol=7)
# for(t in transDates){
#   
#   overallOutput <- rbind(overallOutput,combineBiases(pattern=t))
# }
# finalOutput <- overallOutput[,c(1,2,7,5,6,3,4)]
# colnames(finalOutput) <- c("Site Name","GOES vs. PC","MODIS NDVI vs. PC","MODIS EVI vs. PC",
#                              "MODIS NDVI vs. GOES", "MODIS EVI vs. GOES", "MODIS EVI vs. MODIS NDVI")
# write.table(finalOutput,"GOES_PaperWriting/GOES_PhenologyPaper/PhenologyBias_ALL.csv",sep=",",row.names = FALSE,col.names = TRUE)

#Bias box-and-whisker plot Figure
##I think I need a matrix with one column being the model and one being the bias
#I think that is in each file: 

addBoxPlot <- function(tranName,ylim,label){
  dat <- read.csv(paste("PhenologyBias_GOES_PC_",tranName,".csv",sep=""),header=TRUE)
  dat <- cbind(dat,rep("GOES",nrow(dat)))
  colnames(dat) <- c("siteName","Q1","Q2","Q3","model")
  subDat <- read.csv(paste("PhenologyBias_MODISN_PC_",tranName,".csv",sep=""),header=TRUE)
  subDat <- cbind(subDat,rep("MODIS_NDVI",nrow(subDat)))
  colnames(subDat) <- c("siteName","Q1","Q2","Q3","model")
  subDat
  dat <- rbind(dat,subDat)
  subDat <- read.csv(paste("PhenologyBias_MODISE_PC_",tranName,".csv",sep=""),header=TRUE)
  subDat <- cbind(subDat,rep("MODIS_EVI",nrow(subDat)))
  colnames(subDat) <- c("siteName","Q1","Q2","Q3","model")
  dat <- rbind(dat,subDat)
  
  #dat
  #plot(Q2~model,data=dat,col=GOES.col,pch=20,ylim=ylim,frame.plot=FALSE,xaxt="n",yaxt="n")
  #axis(2,at=c(-30,-20,-10,0,10,20))
  #boxplot(Q2~model, data=dat[dat$model=="MODIS_NDVI",],add=TRUE,col=MODIS.N.col,frame.plot=FALSE,xaxt="n",yaxt="n")
  #boxplot(Q2~model, data=dat[dat$model=="MODIS_EVI",],add=TRUE,col=MODIS.E.col,frame.plot=FALSE,xaxt="n",yaxt="n")
  #axisVals <- c(-20,0,20)
  #axisVals <- seq(ylim[1]+10,ylim[2],20)
  #print(axisVals)
  #axis(2,at=axisVals,labels = label)
  #abline(h=0,col="red",lwd=3)
  returnVal <- dat[dat$model=="MODIS_EVI",][order(abs(dat[dat$model=="MODIS_EVI",]$Q2)),]
  returnV <- cbind(as.character(returnVal[,1]),seq(1,nrow(returnVal),1))
  colnames(returnV) <- c("siteName","ranking")
  print(returnV[order(returnV[,1]),])
  return(returnV[order(returnV[,1]),])
}
#jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_biasPlots_RAWUpdated.jpeg",width=6,height=4.5,units = "in",res=1000)  

par(mfrow=c(2,3),mai=c(0.01,0.4,0.01,0.1))
out <- addBoxPlot(tranName="springTran1",ylim=c(-30,20),label=TRUE)[,c(1,2)]
out <- cbind(out,addBoxPlot(tranName="spring50",ylim=c(-30,20),label=FALSE)[,2])
out <- cbind(out,addBoxPlot(tranName="springTran2",ylim=c(-30,20),label=FALSE)[,2])
out <- cbind(out,addBoxPlot(tranName="fallTran1",ylim=c(-70,80),label=TRUE)[,2])
out <- cbind(out,addBoxPlot(tranName="fall50",ylim=c(-70,80),label=FALSE)[,2])
out <- cbind(out,addBoxPlot(tranName="fallTran2",ylim=c(-70,80),label=FALSE)[,2])
colnames(out) <- c("siteName","springTran1","spring50","springTran2",
                   "fallTran1","fall50","fallTran2")
#write.table(out,file="GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_biasRanking.csv",sep=",",row.names = FALSE,col.names = TRUE)



# addBoxPlot <- function(tranName,ylim,label){
#   dat <- read.csv(paste("PhenologyBias_GOES_PC_",tranName,".csv",sep=""),header=TRUE)[,c(1,3)]
#   #dat <- cbind(dat,rep("GOES",nrow(dat)))
#   colnames(dat) <- c("siteName","Q2")
#   newDat <- dat[order(abs(dat$Q2)),]
# 
#   outVals <-cbind(as.character(newDat[,1]),seq(1,nrow(newDat),1))
#   colnames(outVals) <- c("siteName","GOES")
#   outVals <- outVals[order(outVals[,1]),]
# 
#   dat <- read.csv(paste("PhenologyBias_MODISN_PC_",tranName,".csv",sep=""),header=TRUE)
#   
#   colnames(dat) <- c("siteName","Q2")
#   newDat <- dat[order(abs(dat$Q2)),]
#   
#   outValsNew <-cbind(as.character(newDat[,1]),seq(1,nrow(newDat),1))
#   colnames(outValsNew) <- c("siteName","MODIS_NDVI")
#   outValsNew <- outValsNew[order(outValsNew[,1]),]
#   outVals <- cbind(outVals,outValsNew[,2])
#   
#   dat <- read.csv(paste("PhenologyBias_MODISE_PC_",tranName,".csv",sep=""),header=TRUE)
#   
#   colnames(dat) <- c("siteName","Q2")
#   newDat <- dat[order(abs(dat$Q2)),]
#   
#   outValsNew <-cbind(as.character(newDat[,1]),seq(1,nrow(newDat),1))
#   colnames(outValsNew) <- c("siteName","MODIS_EVI")
#   outValsNew <- outValsNew[order(outValsNew[,1]),]
#   outVals <- cbind(outVals,outValsNew[,2])
#   colnames(outVals) <- c("siteName","GOES","MODIS_NDVI","MODIS_EVI")
#   return(outVals)
#   
# 
# }
# 
# out <- addBoxPlot(tranName="springTran1",ylim=c(-30,20),label=TRUE)
# out <- cbind(out,addBoxPlot(tranName="spring50",ylim=c(-30,20),label=FALSE))
# out <- cbind(out,addBoxPlot(tranName="springTran2",ylim=c(-30,20),label=FALSE))
# out <- cbind(out,addBoxPlot(tranName="fallTran2",ylim=c(-70,80),label=TRUE))
# out <- cbind(out,addBoxPlot(tranName="fall50",ylim=c(-70,80),label=FALSE))
# out <- cbind(out,addBoxPlot(tranName="fallTran1",ylim=c(-70,80),label=FALSE))
# #colnames(out) <- c("siteName","springTran1","spring50","springTran2",
#                    "fallTran1","fall50","fallTran2")
# write.table(out,file="GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_biasRanking.csv",sep=",",row.names = FALSE,col.names = TRUE)
# #dev.off()


######Figure GOES vs MODIS autumn figure


###Examine the rescaled phenology latent states for spring between the different sources
##Will need to read in the data again and calculate the yseq based off of the parameters and xseq


###############Shrubland figure#######################################################
iseq <- c(7,12,13,14)
library("grDevices")
startDay <- 110
endDay <- 455
startDate <- as.Date(startDay,origin="2016-12-31")
endDate <- as.Date(endDay,origin="2016-12-31")
xseq <- seq(startDay,endDay,1)

jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_CompareSourcesShrubland_RAW.jpeg",width=9,height=9,units = "in",res=1000)  

#par(mfrow=c(4,1),mai=c(0.35,1,0.15,0.2))#,mai=c(0.4,0.4,0.2,0.2))
par(mfcol=c(5,2),mai=c(0.30,0.4,0.01,0.2))#,mai=c(0.4,0.4,0.2,0.2))
##Plot Overall GOES
#siteName <- "missouriozarks"
#siteName <- "russellSage"
#siteName <- "shenandoah"
sites <- as.character(read.csv("GOES_Paper_sites.csv",header=TRUE)[iseq,1])
URLs <- as.character(read.csv("GOES_Paper_sites.csv",header=TRUE)[iseq,4])
GOES.ks <- matrix(nrow=0,ncol=3)
MODISN.ks <- matrix(nrow=0,ncol=3)
MODISE.ks <- matrix(nrow=0,ncol=3)
PC.ks <- matrix(nrow=0,ncol=3)

for(i in 1:length(sites)){
  siteName <- sites[i]
  print(siteName)
  BasicInFileName <- paste(siteName,"_GOES_noon_varBurn.RData",sep="")
  noonDat <- read.csv(paste("GOES_NDVI_",siteName,"_",startDate,"_",endDate,"_noon.csv",sep=""),header=FALSE)
  #noonDat <- read.csv("GOES_NDVI_HarvardForest_2017-07-01_2018-06-30_noon.csv",header=FALSE)
  for(l in 1:length(noonDat)){
    if(noonDat[1,l]<startDay){
      noonDat[1,l] <- noonDat[1,l]+365
    }
  }
  noon.dates <- as.Date(as.numeric(noonDat[1,]),origin="2016-12-31")
  
  load(BasicInFileName)
  var.mat<-data.frame(as.matrix(GOES.md.out))
  colnames(var.mat)
  
  GOES.k <- var.mat$k
  
  ci.Overall <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale = FALSE)
  ci.Overall2 <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale = TRUE)
  
  
  plot(noon.dates,noonDat[2,],pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0.1,0.9))
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col=GOES.col) #"#a6cee3"
  
  points(noon.dates,noonDat[2,],pch=20,col="black")
  
  print("Done with GOES")
  ##Plot PhenoCam
  data.PC = PC_data(siteName=siteName,URL=URLs[i],startDate=startDate,endDate=endDate)
  inputFileName <- paste(siteName,"_110_455_PC_varBurn.RData",sep="")
  inputFileName <- paste(siteName,"_PC_varBurn.RData",sep="")
  #load("luckyHills_110_455_PC_varBurn.RData")
  load(inputFileName)
  
  var.mat<-data.frame(as.matrix(PC.md.out))
  PC.k <- var.mat$k
  
  
  ci.PC <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale=FALSE)
  ci.PC2 <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale=TRUE)
  PC.x <- data.PC$x#[data.PC$x>(startDay-1)]
  PC.y <- data.PC$y#[data.PC$x>(startDay-1)]
  PC.dates <- as.Date(PC.x,origin="2016-12-31")
  
  
  plot(PC.dates,PC.y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0.3,0.46))
  #ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="#a6cee3")
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.PC[1,],ci.PC[3,],col=PC.col)#"#1f78b4")
  points(PC.dates,PC.y,pch=20)
  print("Done with PC")
  ##MODIS NDVI
  data.MODIS.N = MODIS_data(siteName=siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="NDVI")
  inputFileName <- paste(siteName,"_110_455_MODIS_DQF_NDVI_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-data.frame(as.matrix(MODIS.N.md.out))
  MODISN.k <- var.mat$k
  
  ci.MODIS.N <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale = FALSE)
  ci.MODIS.N2 <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale = TRUE)
  
  MODIS.dates <- as.Date(data.MODIS.N$x,origin="2016-12-31")
  #ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="#a6cee3")
  plot(MODIS.dates,data.MODIS.N$y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0,1))
  
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.MODIS.N[1,],ci.MODIS.N[3,],col=MODIS.N.col)
  points(MODIS.dates,data.MODIS.N$y,pch=20)
  print("Done with MODIS NDVI")
  ##MODIS EVI
  data.MODIS.E = MODIS_data(siteName=siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="EVI")
  inputFileName <- paste(siteName,"_110_455_MODIS_DQF_EVI_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-data.frame(as.matrix(MODIS.E.md.out))
  MODISE.k <- var.mat$k
  
  ci.MODIS.E <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale = FALSE)
  ci.MODIS.E2 <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale = TRUE)
  
  MODIS.dates <- as.Date(data.MODIS.E$x,origin="2016-12-31")
  #ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="#a6cee3")
  plot(MODIS.dates,data.MODIS.E$y,pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n',ylim=c(0,0.8))
  
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.MODIS.E[1,],ci.MODIS.E[3,],col=MODIS.E.col)
  points(MODIS.dates,data.MODIS.E$y,pch=20)
  print("Done with MODIS EVI")
  
  plot(dat.dates,rep(NA,length(dat.dates)),ylim=c(-0.15,1),pch=20,col="black",ylab="",cex.axis=2,cex.lab=2,xlab="",main="",bty = 'n')
  
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.MODIS.N2[1,],ci.MODIS.N2[3,],col=MODIS.N.col)
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.MODIS.E2[1,],ci.MODIS.E2[3,],col=adjustcolor(MODIS.E.col,0.6))
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall2[1,],ci.Overall2[3,],col=adjustcolor(GOES.col,0.6))
  ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.PC2[1,],ci.PC2[3,],col=adjustcolor(PC.col,0.6))#"#1f78b4")
  
  # polygon(x=c(as.Date(quantile(GOES.Tran,0.025),origin="2016-12-31"),
  #             as.Date(quantile(GOES.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(GOES.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(GOES.Tran,0.025),origin="2016-12-31")),
  #         y=c(-0.1,-0.1,-0.02,-0.02),col=adjustcolor(GOES.col,0.5),border=NA)
  # lines(x=c(as.Date(quantile(GOES.Tran,0.5),origin="2016-12-31"),as.Date(quantile(GOES.Tran,0.5),origin="2016-12-31"))
  #       ,y=c(-0.1,-0.02),col=GOES.col,lwd=3)
  # 
  # polygon(x=c(as.Date(quantile(PC.Tran,0.025),origin="2016-12-31"),
  #             as.Date(quantile(PC.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(PC.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(PC.Tran,0.025),origin="2016-12-31")),
  #         y=c(-0.1,-0.1,-0.02,-0.02),col=adjustcolor(PC.col,0.5),border=NA)
  # lines(x=c(as.Date(quantile(PC.Tran,0.5),origin="2016-12-31"),as.Date(quantile(PC.Tran,0.5),origin="2016-12-31"))
  #       ,y=c(-0.1,-0.02),col=PC.col,lwd=3)
  # 
  # polygon(x=c(as.Date(quantile(MODISN.Tran,0.025),origin="2016-12-31"),
  #             as.Date(quantile(MODISN.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(MODISN.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(MODISN.Tran,0.025),origin="2016-12-31")),
  #         y=c(-0.1,-0.1,-0.02,-0.02),col=adjustcolor(MODIS.N.col,0.5),border=NA)
  # lines(x=c(as.Date(quantile(MODISN.Tran,0.5),origin="2016-12-31"),as.Date(quantile(MODISN.Tran,0.5),origin="2016-12-31"))
  #       ,y=c(-0.1,-0.02),col=MODIS.N.col,lwd=3)
  # 
  # 
  # polygon(x=c(as.Date(quantile(MODISE.Tran,0.025),origin="2016-12-31"),
  #             as.Date(quantile(MODISE.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(MODISE.Tran,0.975),origin="2016-12-31"),
  #             as.Date(quantile(MODISE.Tran,0.025),origin="2016-12-31")),
  #         y=c(-0.1,-0.1,-0.02,-0.02),col=adjustcolor(MODIS.E.col,0.5),border=NA)
  # lines(x=c(as.Date(quantile(MODISE.Tran,0.5),origin="2016-12-31"),as.Date(quantile(MODISE.Tran,0.5),origin="2016-12-31"))
  #       ,y=c(-0.1,-0.02),col=MODIS.E.col,lwd=3)
  GOES.ks <- rbind(GOES.ks,quantile(GOES.k,c(0.025,0.5,0.975)))
  MODISN.ks <- rbind(MODISN.ks,quantile(MODISN.k,c(0.025,0.5,0.975)))
  MODISE.ks <- rbind(MODISE.ks,quantile(MODISE.k,c(0.025,0.5,0.975)))
  PC.ks <- rbind(PC.ks,quantile(PC.k,c(0.025,0.5,0.975)))
  
  
}



dev.off()



par(pty="s")
plot(PC.CIs.TranSs1,GOES.CIs.TranSs1,pch=20,ylim=c(0,40),xlim=c(0,40),cex=2)
points(PC.CIs.TranSs1,MODIS.E.CIs.TranSs1,pch=20,ylim=c(0,40),xlim=c(0,40),col=MODIS.E.col,cex=2)
points(PC.CIs.TranSs1,MODIS.N.CIs.TranSs1,pch=20,ylim=c(0,40),xlim=c(0,40),col=MODIS.N.col,cex=2)
abline(0,1,col="red")
