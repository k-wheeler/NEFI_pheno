##Script to classify the pattern for each diurnal day
##Probably need to write functions for each classification, loop over all of the days and input each day's diurnal data 
##into each of the functions, which will return true or false depending on if the classification holds or not

library("rjags")
library("suncalc")
days <- seq(as.Date("2017-07-01"),as.Date("2018-06-30"),"day")
#iseq <- as.character(c(seq(1,333,1),seq(348,364,1)))
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
s <- 1

for(s in 1:nrow(siteData)){
  PFT <- as.character(siteData$PFT[s])
  if(PFT=="DB"){
    siteName <- as.character(siteData$siteName[s])
    TZ <- as.character(siteData$TZ[s])
    print(siteName)
    dayObs <- numeric()
    clfts <- character() #List of diurnal pattern classifications
    lat <- as.numeric(siteData$Lat[s])
    long <- as.numeric(siteData$Long[s])
    for(d in 1:length(days)){
      ##Load in the data
      dy <- as.character(format(days[d],"%j"))
      year <- lubridate::year(days[d])
      fileName <- paste("dailyNDVI_GOES/",siteName,"_GOES_diurnal_",year,dy,".csv",sep="")
      print(fileName)
      if(file.exists(fileName)){
        dat <- read.csv(fileName,header=FALSE)
        dayObs <- c(dayObs,dy)
        data <- list()
        print(dim(dat))
        data$x <- as.numeric(dat[3,])
        data$y <- as.numeric(dat[2,])
        totalLength <- c(totalLength,length(na.omit(data$y)))
        ##Calculate noise:
        #plot(data$x,data$y,pch=20)
        midWindow <- data$y[data$x>10 &data$x<14]
        winCI <- as.numeric(quantile(midWindow,0.975,na.rm=TRUE)-quantile(midWindow,0.025,na.rm=TRUE))
        #points(data$x[data$x>10 &data$x<14],midWindow,col="cyan",pch=20)
        
        winLength <-length(na.omit(midWindow))
        
        if(TZ==6){
          solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(dy),origin=as.Date(paste((year-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/Chicago"))$solarNoon
          solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
        }
        else if(TZ==5){
          solarNoon <- (getSunlightTimes(date=as.Date(as.numeric(dy),origin=as.Date(paste((year-1),"-12-31",sep=""))),lat=lat,lon=long,keep="solarNoon",tz="America/New_York"))$solarNoon
          solarNoonTime <- lubridate::hour(solarNoon)+(lubridate::minute(solarNoon)/60)
        }
        
        fitFileName <- paste(siteName,"_",dy,"_varBurn6.RData",sep="")
        if(file.exists(fitFileName)){
          load(fitFileName)
          out.mat <- data.frame(as.matrix(var.burn))
          c <- out.mat$c
          rndNums <- sample(1:length(c),10000,replace=T)
          c <- c[rndNums]
          fitCI <- as.numeric(quantile(c,0.975))-as.numeric(quantile(c,0.025))
        }
        else{
          fitCI <- NA
          clf <- "NoFit"
        }
        ##Classify day:
        if(totalLength<10){
          clf <- "NoData"
        }
        else{
          if(winLength<5){
            if(fitCI<0.1){
              clf <- "NoWindow_tightFit"
            }
            else{
              clf <- "NoWindow_wideFit"
            }
          }
          else{
            if(fitCI<0.1){
              if(winCI<0.1){
                clf <- "LowNoise_tightFit"
              }
              else{
                clf <- "HighNoise_tightFit"
              }
            }
            else{
              if(winCI<0.1){
                clf <- "LowNoise_wideFit"
              }
              else{
                clf <- "HighNoise_wideFit"
              }
            }
          }
        } ##Ends the totalLength>=10
        
      } ##Ends if file.exists()
      else{
        clf <- "NoData"
      }
      clfts <- c(clfts,clf)
      dayObs <- c(dayObs,dy)
    }
    #output <- cbind(dayObs,totalLength,winLength,winCI,fitCI)
    output <- cbind(dayObs,clfts)
    outFileName <- paste(siteName,"_diurnal_classifications.csv",sep="")
    write.table(output,file=outFileName,sep=",",row.names=FALSE,col.names = TRUE)
    
  }
}
#   
#   
# dat <- read.csv(outFileName,header=TRUE) 
# 
# nrow(dat[dat$winLength==0,])
# nrow(dat[dat$totalLength>10,])
# newDat <- dat[dat$totalLength>10,]
# 
# nrow(newDat[newDat$fitCI<0.1,])
# nrow(newDat[newDat$winCI<0.1,])
# hist(newDat$winCI)
# hist(newDat$winLength[newDat$winLength!=0])
# newDat2 <- newDat[newDat$winLength>5,]
# nrow(newDat2[newDat2$winCI<0.05,])
# 
# nrow(newDat[newDat$fitCI<0.05,])
# 
# avg.sd <- sqrt(var(as.numeric(midday[2,]),na.rm=TRUE))
# avg.bot <- avg-1.96*(avg.sd/sqrt(length(na.omit(as.numeric(midday[2,])))))
# 
# ##A) Number with >10 observations
# 
# ##B) From A Number with >5 observations within the mid-window range
# 
# ##c) From B, number with winCI < 0.05 <- low noise
# 
# ##
# 
# ##Classifications
# missingMid_tightFit <- newDat[(newDat$winLength<5&newDat$fitCI<0.1),]
# 
# highNoise_tightFit <- newDat[(newDat$winCI>0.1 & newDat$fitCI<0.1 & !is.na(newDat$winCI)),]
# 
# highNoise_loseFit <- newDat[(newDat$winCI>0.1 & newDat$fitCI>0.1 & !is.na(newDat$winCI)),]
# 
# lowNoise_tightFit <- newDat[(newDat$winCI<0.1 & newDat$fitCI<0.1 & !is.na(newDat$winCI)),]
# 
# lowNoise_loseFit <- newDat[(newDat$winCI<0.1 & newDat$fitCI>0.1 & !is.na(newDat$winCI)),]
# 
# 



  
