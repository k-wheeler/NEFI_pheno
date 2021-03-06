#!/usr/bin/env Rscript

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
library("ncdf4")
library(plyr)
library("PhenologyBayesModeling")
library(doParallel)
library("rjags")
library("runjags")
n.cores <- 5

#register the cores.
#registerDoParallel(cores=n.cores)

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
iseq <- c(seq(1,6),8,9,10,seq(15,20))
iseq <- seq(1,6)
#output <- foreach(i = iseq) %dopar% {
for(i in iseq){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  PFT <- as.character(siteData$PFT[i])
  print(PFT)
  if(PFT=="DB"){
    #diurnalFits <- dir(path="diurnalFits",pattern=siteName)
    diurnalFits <- intersect(dir(pattern="varBurn6.RData"),dir(pattern=siteName))
    c.vals <- numeric()
    Q1.vals <- numeric()
    Q2.vals <- numeric()
    Q3.vals <- numeric()
    days <- numeric()
    counts <- numeric()
    outDataFile <- paste(siteName,"_diurnal6FitData.RData",sep="")
    if(!file.exists(outDataFile)){
      for(i in 1:length(diurnalFits)){
        print(diurnalFits[i])
        #load(paste("diurnalFits/",diurnalFits[i],sep=""))
        load(diurnalFits[i])
        if(typeof(var.burn)!=typeof(FALSE)){
          out.mat <- as.matrix(var.burn)
          print(colnames(out.mat))
          c <- mean(out.mat[,2])
          Q1 <- as.numeric(quantile(out.mat[,2],0.025))
          Q2 <- as.numeric(quantile(out.mat[,2],0.5))
          Q3 <- as.numeric(quantile(out.mat[,2],0.975))
          #prec <- as.numeric(quantile(out.mat[,2],0.975))-as.numeric(quantile(out.mat[,2],0.025))
          dy <- strsplit(diurnalFits[i],"_")[[1]][2]
          if(dy<182){
            yr <- 2018
          }
          else{
            yr <- 2017
          }
          print(paste("_GOES_diurnal_",dy,".csv",sep=""))
          dayDataFile <- intersect(dir(path="dailyNDVI_GOES",pattern=paste("_GOES_diurnal_",yr,dy,".csv",sep="")),dir(path="dailyNDVI_GOES",pattern=siteName))
          print(dayDataFile)
          dayData <- read.csv(paste("dailyNDVI_GOES/",dayDataFile,sep=""),header=FALSE)
          ct <- length(dayData[2,][!is.na(dayData[2,])])
          
          c.vals <- c(c.vals,c)
          Q1.vals <- c(Q1.vals,Q1)
          Q2.vals <- c(Q2.vals,Q2)
          Q3.vals <- c(Q3.vals,Q3)
          #prec.vals <- c(prec.vals,prec)
          counts <- c(counts,ct)
          days <- c(days,dy)
          
        }
      }
      data <- list()
      for(i in 1:length(days)){
        if(days[i]<182){
          days[i] <- as.numeric(days[i]) + 365
        }
      }
      data$x <- as.numeric(days)
      data$y <- as.numeric(c.vals)
      data$Q1<- as.numeric(Q1.vals)
      data$Q2<- as.numeric(Q2.vals)
      data$Q3<- as.numeric(Q3.vals)
      data$n <- length(data$x)
      data$size <- as.numeric(counts)
      print(dim(data$x))
      print(dim(data$y))
      print(data$x)
      save(data,file=outDataFile)
      print(paste("Done with creating Data for", siteName))
    }
  }
}
