#!/usr/bin/env Rscript

#install.packages("devtools")
#library("devtools")
#install_github("EcoForecast/ecoforecastR")
library("ecoforecastR")
library("rjags")
library("suncalc")
library("doParallel")

n.cores <- 5

#register the cores.
registerDoParallel(cores=n.cores)

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

# calSatAlt <- function(lat,long){
#   ##Define constants
#   satLat <- 0
#   satLong <- -75.2
#   
#   R <- 6378140 #(in m from class notes)
#   Rheight <- 35800*1000 #https://noaasis.noaa.gov/NOAASIS/ml/genlsatl.html
#   Rsat <- R + Rheight
#   
#   term1 <- Rsat*(cos(satLong-long)*cos(satLat)*cos(lat)+sin(satLat)*sin(lat))-R
#   term2 <- Rsat**2+R**2-2*R*Rsat*(cos(satLong - long)*cos(satLat)*cos(lat)+sin(satLat)*sin(lat))
#   product <- term1 * term2**(-1/2)
#   return(1.5708-acos(product))
# }
# calSatAzm <- function(lat,long){
#   ##Define constants
#   satLat <- 0
#   satLong <- -75.2
#   
#   R <- 6378140 #(in m from class notes)
#   Rheight <- 35800*1000 #https://noaasis.noaa.gov/NOAASIS/ml/genlsatl.html
#   Rsat <- R + Rheight
#   
#   num <- sin(satLong-long)*cos(satLat)
#   den <- cos(satLong - long)*cos(satLat)*sin(lat)-sin(satLat)*cos(lat)
#   quot <- num/den
#   return(atan(quot))
# }

#outputFileName <- "ALL_DiurnalFits.pdf"
#pdf(file=outputFileName,width=45,height=40)
xseq <- seq(0,25,0.1)
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
#iseq <- c(seq(1,6),seq(8,11),seq(15,20))
#iseq <- c(seq(4,6),seq(8,11),seq(15,20))
#iseq <- c(8,9)
#sseq <- c(seq(1,6),8,9,10,11,seq(15,20))
classData <- read.csv("HarvardForest_diurnal_characteristics.csv",header=TRUE)
s <- 1

#output <- foreach(s = sseq) %dopar% {
#for(s in iseq){
  siteName <- as.character(siteData[s,1])
  outputFileName <- paste(siteName,"_ALL_DiurnalFits6.pdf",sep="")
  lat <- as.numeric(siteData$Lat[s])
  long <- as.numeric(siteData$Long[s])
  pdf(file=outputFileName,width=20,height=20)
  par(mfrow=c(5,5))
  diurnalFiles <- intersect(dir(pattern="varBurn6.RData"),dir(pattern=siteName))
  #sunAngleFile <- paste(siteName,"_sunAngles.csv",sep="")
  #sunAngles <- read.csv(sunAngleFile,header=TRUE,sep=" ")
  
  #satAlt <- calSatAlt(lat=lat,long=long)
  #satAzm <- calSatAzm(lat=lat,long=long)
  # timeStr <- character()
  # for(a in 1:nrow(sunAngles)){
  #   if(a %% 1000 == 0){
  #     print(a)
  #     print(paste(yr,dy,hr,mn,sep=""))
  #   }
  #   yr <- as.character(lubridate::year(sunAngles[a,1]))
  #   dy <- as.character(format(as.Date(sunAngles[a,1]),"%j"))
  #   hr <- as.character(lubridate::hour(sunAngles[a,1]))
  #   mn <- as.character(lubridate::minute(sunAngles[a,1]))
  #   if(as.numeric(dy)<10){
  #     dy <- paste("00",as.numeric(dy),sep="")
  #   }
  #   else if(as.numeric(dy)<100){
  #     dy <- paste("0",as.numeric(dy),sep="")
  #   }
  #   if(as.numeric(hr)<10){
  #     hr <- paste("0",as.numeric(hr),sep="")
  #   }
  #   if(as.numeric(mn)<10){
  #     mn <- paste("0",as.numeric(mn),sep="")
  #   }
  # 
  #   timeStr <- c(timeStr,paste(yr,dy,hr,mn,sep=""))
  # }
  # sunAngles <- cbind(sunAngles,timeStr)
  
  for(i in 1:length(diurnalFiles)){
    load(diurnalFiles[i])
    dy <- strsplit(diurnalFiles[i],"_")[[1]][2]
    print(diurnalFiles[i])
    print(dy)
    print(as.numeric(dy))
    if(as.numeric(dy)<182){
      yr <- "2018"
    }
    if(as.numeric(dy)>181){
    #else{
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
    ##Check for hot spot potentials:
    # for(c in 1:ncol(dat)){
    #   which(as.character(sunAngles$timeStr)==as.character(dat[1,c]))
    # }
    # as.character(dat[1,c])
    
    if(typeof(var.burn)==typeof(FALSE)){
      print(paste(diurnalFiles[i], " did not converge",sep=""))
      #plot(as.numeric(dat[3,]),as.numeric(dat[2,]),main=paste("Didn't",diurnalFiles[i],sep=" "),xlab="Time",ylab="NDVI",ylim=c(0,1),xlim=c(0,25))
    }
    else{
      out.mat <- data.frame(as.matrix(var.burn))
      a <- out.mat$a
      rndNums <- sample(1:length(a),10000,replace=T)
      a <- a[rndNums]
      c <- out.mat$c[rndNums]
      #alp <- out.mat$alp[rndNums]
      #bet <- out.mat$bet[rndNums]
      #p.cloud <- out.mat$p.cloud[rndNums]
      k <- out.mat$k[rndNums]
      solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(i),origin=as.Date(paste((as.numeric(yr)-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/Chicago"))$solarNoon
      solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
      ycred <- matrix(0,nrow=10000,ncol=length(xseq))
      for(g in 1:10000){
        Ey <- diurnalExp(a=a[g],c=c[g],k=k[g],xseq=xseq)
        ycred[g,] <- Ey
      }
      ci <- apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE)
      #plot(x=list(),y=list(),main=diurnalFiles[i],xlab="Time",ylab="NDVI",ylim=c(0,1),xlim=c(0,25))
      if(length(na.omit(as.numeric(dat[2,])))>10){
        ##Calculate 
        classDat <- classData[classData$dayObs==dy,]
        print(classDat)
        if(classDat$winLength<10){
          if(classDat$fitCI<0.1){
            clf <- "No window, tight fit"
          }
          else{
            clf <- "No window, wide fit"
          }
        }
        else{
          if(classDat$fitCI<0.1){
            if(classDat$winCI<0.1){
              clf <- "Low noise, tight fit"
            }
            else{
              clf <- "High noise, tight fit"
            }
          }
          else{
            if(classDat$winCI<0.1){
              clf <- "Low noise, wide fit"
            }
            else{
              clf <- "High noise, wide fit"
            }
          }
        }
        
        plot(as.numeric(dat[3,]),as.numeric(dat[2,]),main=paste(siteName,dy,clf),xlim=c(0,25),pch=20,cex=2,ylim=c(0,1.2))
        ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
        lines(xseq,ci[2,],col="black")
        points(as.numeric(dat[3,]),as.numeric(dat[2,]),pch=20,cex=2)
        #abline(v=12,col="red")
        abline(v=solarNoonTime,col="purple")
      #}
    }
  #}
    }
  }
  dev.off()
#}
    
#plot(density(rbeta(10000,1,25)))
#abline(v=0.13,col="red")
    
#sum(rbeta(1000000,1,17)<0.13)/10000
