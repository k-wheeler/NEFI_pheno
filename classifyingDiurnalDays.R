##Script to classify the pattern for each diurnal day
##Probably need to write functions for each classification, loop over all of the days and input each day's diurnal data 
  ##into each of the functions, which will return true or false depending on if the classification holds or not

library("rjags")
library("suncalc")
days <- seq(as.Date("2017-07-01"),as.Date("2018-06-30"),"day")
#iseq <- as.character(c(seq(1,333,1),seq(348,364,1)))
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
s <- 1
dayObs <- numeric()
winCI <- numeric()
winLength <- numeric()
totalLength <- numeric()
fitCI <- numeric()
#for(s in iseq){
  siteName <- as.character(siteData$siteName[s])
  TZ <- as.character(siteData$TZ[s])
  lat <- as.numeric(siteData$Lat[s])
  long <- as.numeric(siteData$Long[s])
  print(siteName)
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
      winCI <- c(winCI,as.numeric(quantile(midWindow,0.975,na.rm=TRUE)-quantile(midWindow,0.025,na.rm=TRUE)))
      #points(data$x[data$x>10 &data$x<14],midWindow,col="cyan",pch=20)
      
      winLength <- c(winLength,length(na.omit(midWindow)))
      
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
        fitCI <- c(fitCI,as.numeric(quantile(c,0.975))-as.numeric(quantile(c,0.025)))
      }
      else{
        fitCI <- c(fitCI,NA)
      }
    }
  }
  output <- cbind(dayObs,totalLength,winLength,winCI,fitCI)
  outFileName <- paste(siteName,"_diurnal_characteristics.csv",sep="")
  write.table(output,file=outFileName,sep=",",row.names=FALSE,col.names = TRUE)
  
#}