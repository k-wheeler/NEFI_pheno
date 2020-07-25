s <- 1
days <- seq(1,181)
#' Creates the diurnal jags model object
#'
#' @param siteName The site name used for file naming
#' @param data List that includes the data to be fitted (e.g. data$x and data$y) where x is the vector of hours of the day (minutes and seconds included in the decimal) y is the vector of NDVI values
#' @export
#' @import rjags
#' @import runjags
#' @import coda
createDiurnalModel2 <- function(siteName,data){
  print("entered model")
  nchain <-  5
  data$alpha.c <- 1#2
  data$beta.c <- 1#1.5
  print(data$beta.c)
  data$s1 <- 0.001
  data$s2 <- 0.00001
  
  data$mean.a <- 0.0009
  data$p.a <- 0.0003
  data$mean.k <- 12
  data$p.k <- 1/(1**2)
  data$n <- length(data$x)
  
  inits <- list()
  middleSection <- numeric()
  for(i in 1:length(data$x)){
    if(data$x[i]<15 && data$x[i]>10){
      middleSection <- c(middleSection,data$y[i])
    }
  }
  #print(length(middleSection))
  if(length(na.omit(middleSection))<1){
    middleSection <- 0.50
  }
  print(mean(middleSection,na.rm=TRUE))
  for(n in 1:nchain){
    c <- rnorm(1,mean(middleSection,na.rm = TRUE),0.1)
    while(c<0.1 || c>mean(middleSection,na.rm=TRUE)+1.5){
      c <- rnorm(1,mean(middleSection,na.rm=TRUE),0.05)
    }
    inits[[n]] <- list(c=c,a=rnorm(1,0.001,0.0001),k=rnorm(1,14.5,0.5))
  }
  print(inits)
  print("finished defining data")
  
  DB_model_MM <- "
  model{
  ##priors
  #TranL ~ dnorm(mean.TranL,p.Tran) ##S for spring
  #bL ~ dnorm(mean.bL,p.b)
  #TranR ~ dnorm(mean.TranR,p.Tran)  ##F for fall/autumn
  #bR ~ dnorm(mean.bR,p.b)
  a ~ dnorm(mean.a,p.a) I(0,)
  c ~ dbeta(alpha.c,beta.c)
  k ~ dnorm(mean.k,p.k)
  prec ~ dgamma(s1,s2)
  alp ~ dunif(1,100)
  bet ~ dunif(1,100)
  p.cloud ~ dunif(0,1)
  
  for(i in 1:n){
  muL[i] <- -a * exp(-1 * (x[i]-k)) + c + a
  muR[i] <- -a * exp((x[i]-k)) + c + a
  
  f[i] <- ifelse(x[i]>k,muR[i],muL[i])   #change point process model
  
  y[i] ~ dnorm(mu[i],prec)   ##data model
  is.cloudy[i] ~ dbern(p.cloud)
  trans[i] ~ dbeta(alp,bet)
  mu[i] <- is.cloudy[i] * trans[i]*f[i] + (1-is.cloudy[i]) * f[i]
  
  }
  }
  "
  
  j.model   <- jags.model(file = textConnection(DB_model_MM),
                          inits = inits,
                          data = data,
                          n.chains=nchain)
  return(j.model)
  
}


diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  bk <- which(round(xseq,digits=1)==k)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  return(c(left,right))
}

library(suncalc)
library(GOESDiurnalNDVI)
library("ecoforecastR")
library("rjags")
library("runjags")
library("doParallel")

n.cores <- 15

#register the cores.
registerDoParallel(cores=n.cores)

siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
savePath <- paste(getwd(),"/PhenologyForecastData/GOES_DiurnalFits/",sep="")

siteName <- as.character(siteData$siteName[s])
lat <- as.numeric(siteData$Lat[s])
long <- as.numeric(siteData$Long[s])
TZ <- as.character(siteData$TZ[s])
if(TZ==5){
  TZ_name <- "America/New_York"
}else if(TZ==6){
  TZ_name <- "America/Chicago"
}else if(TZ_name==""){
  print("TZ_name unknown. Please enter it")
  return(FALSE)
}
PFT <- as.character(siteData$PFT[s])
print(siteName)    

#days <- seq(1,365)
year <- 2018

for(i in 1:length(days)){
  if(as.numeric(days[i])<10){
    days[i] <- paste("00",as.character(days[i]),sep="")
  }
  else if(as.numeric(days[i])<100){
    days[i] <- paste("0",days[i],sep="")
  }
}

output <- foreach(day = days) %dopar% {
#for(day in days){
  #fileName <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="") #For 2019 and second half of 2018
  fileName <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/",siteName,"_GOES_diurnal_",year,day,".csv",sep="")
  print(fileName)

  dayData <- read.csv(fileName,header=FALSE)
  if(ncol(dayData)>1){ 
    date.val <- as.Date(as.numeric(day),origin=as.Date(paste(as.character(as.numeric(year)-1),"-12-31",sep="")))
    
    suntimes <- getSunlightTimes(date=date.val,lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),tz = TZ_name)
    dawnTime <- lubridate::hour(suntimes$nauticalDawn)+(lubridate::minute(suntimes$nauticalDawn)/60)
    duskTime <- lubridate::hour(suntimes$nauticalDusk)+(lubridate::minute(suntimes$nauticalDusk)/60)
    
    for(c in 1:ncol(dayData)){
      if(!is.na(dayData[2,c])){
        if(round(dayData[2,c],digits=4)==0.6040){
          dayData[2,c] <- NA
        }else if(dayData[2,c]>0.9){
          dayData[2,c] <- NA
        }else if(dayData[2,c]<0){
          dayData[2,c] <- NA
        }
      }
    }
    
    dayData[2,as.numeric(dayData[3,])<(as.numeric(dawnTime)+1.5)] <- NA
    dayData[2,as.numeric(dayData[3,])>(as.numeric(duskTime)-1.5)] <- NA
    data <- list(x=as.numeric(dayData[3,]),y=as.numeric(dayData[2,]))
    modelFitFileName <- paste(savePath,siteName,"_",year,day,"_varBurnFilter.RData",sep="")
    if(!file.exists(modelFitFileName)){
      j.model <- createDiurnalModel2(siteName=siteName,data=data)
      var.burn <- runMCMC_Model(j.model=j.model,variableNames=c("a","c","k","prec"),
                                baseNum=40000,iterSize =80000,maxGBR = 1.3)

      if(typeof(var.burn)!=typeof(FALSE)){
        out.mat <- as.matrix(var.burn)
        thinAmount <- round(nrow(out.mat)/5000,digits=0)
        var.burn <- window(var.burn,thin=thinAmount)
        save(var.burn,file=modelFitFileName)
      }
    }
  }
}



