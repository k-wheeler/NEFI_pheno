##' Create a Bayes Model for a deciduous broadleaf site using all data sources
##'
##' @param siteName Site Name
##' @param URL PhenoCam network URL
##' @param lat latitude of site in degrees
##' @param long longitude of site in degrees
##' @param startDay the day of year since 2017-01-01 to start the model
##' @param endDay the day of year since 2017-01-01 to end the model
##' @param niter the maximum number of iterations you want to give the model to converge within
##' @param seasonOrder The order of seasons (default is "FS" for fall and then spring)
##' @param dataDirectory The data directory
##' @import rjags
##' @import runjags
##' @export
createDBModel_ALL <- function(siteName="",URL="",niter=100000,startDay,endDay,lat,long,TZ=5,seasonOrder="FS",dataDirectory) {
  nchain = 5
  data <- list()

  p.data <- PC_data(siteName=siteName,URL=URL,startDate=as.Date("2018-07-01"),endDate=as.Date("2019-06-30"))
  mn.data = MODIS_data(siteName=siteName,lat=lat,long=long,startDate = startDate,endDate = endDate,metric="NDVI",lastYear = 2019,dataDirectory=dataDirectory)

  me.data = MODIS_data(siteName=siteName,lat=lat,long=long,startDate = startDate,endDate = endDate,metric="EVI",lastYear=2019,dataDirectory=dataDirectory)


  ##Need to format by adding NAs to missing days
  xseq <- seq(startDay,endDay,1)
  newP <- numeric()
  newME <- numeric()
  newMN <- numeric()
  for(i in 1:length(xseq)){
    if(xseq[i]%in% p.data$x){
      newP <- c(newP,p.data$y[p.data$x==xseq[i]])
    }else{
      newP <- c(newP,NA)
    }
    if(xseq[i]%in% mn.data$x){
      newMN <- c(newMN,mn.data$y[mn.data$x==xseq[i]])
    }else{
      newMN <- c(newMN,NA)
    }
    if(xseq[i]%in% me.data$x){
      newME <- c(newME,me.data$y[me.data$x==xseq[i]])
    }else{
      newME <- c(newME,NA)
    }
  }

  ###
  ##Need to rescale the data
  rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_PC.csv",sep="")
  rescaleData <- read.csv(rescaleFile,header=TRUE)
  cMeans.p <- rescaleData$cMeans.p
  dMeans.p <- rescaleData$dMeans.p
  rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_MN.csv",sep="")
  rescaleData <- read.csv(rescaleFile,header=TRUE)
  cMeans.mn <- rescaleData$cMeans.mn
  dMeans.mn <- rescaleData$dMeans.mn
  rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_ME.csv",sep="")
  rescaleData <- read.csv(rescaleFile,header=TRUE)
  cMeans.me <- rescaleData$cMeans.me
  dMeans.me <- rescaleData$dMeans.me

  data$p <- rescaleObs(times=as.Date(xseq,origin=as.Date("2017-12-31")),vals=newP,partialStart=FALSE,cVals=cMeans.p,dVals=dMeans.p)
  data$mn <- rescaleObs(times=as.Date(xseq,origin=as.Date("2017-12-31")),vals=newMN,partialStart=FALSE,cVals=cMeans.mn,dVals=dMeans.mn)
  data$me <- rescaleObs(times=as.Date(xseq,origin=as.Date("2017-12-31")),vals=newME,partialStart=FALSE,cVals=cMeans.me,dVals=dMeans.me)

  data$n <- length(xseq)
  data$mean.TranS <- 475
  data$p.Tran <- 1/(40**2)
  data$mean.bS <- -0.10
  data$mean.TranF <- 300
  data$mean.bF <- 0.10
  data$p.b <-  1/(0.05**2)
  data$mean.k <- 365
  data$p.k <- 1/(1**2)
  data$s1.PC <- 0.5
  data$s2.PC <- 0.2
  data$s1.MN <- 0.5
  data$s2.MN <- 0.2
  data$s1.ME <- 0.5
  data$s2.ME <- 0.2
  data$a.d <- 1 #minimum
  data$b.d <-10
  data$a.c <- 10
  data$b.c <-1
  data$x <- xseq
  #plot(data$x,data$p,pch=20)
  #points(data$x,data$mn,pch=20,col="blue")
  #points(data$x,data$me,pch=20,col="red")

  if(seasonOrder=="FS"){
    DB_model <- "
    model{
    ##priors
    TranS ~ dnorm(mean.TranS,p.Tran) ##S for spring
    bS ~ dnorm(mean.bS,p.b)
    TranF ~ dnorm(mean.TranF,p.Tran)  ##F for fall/autumn
    bF ~ dnorm(mean.bF,p.b)
    d ~ dbeta(a.d,b.d)
    c ~ dbeta(a.c,b.c)
    k ~ dnorm(mean.k,p.k)
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.MN ~ dgamma(s1.MN,s2.MN)
    p.ME ~ dgamma(s1.ME,s2.ME)

    for(i in 1:n){
    muF[i] <- c/(1+exp(bF*(x[i]-TranF)))+d ##process model for fall
    muS[i] <- c/(1+exp(bS*(x[i]-TranS)))+d ##process model for Spring
    mu[i] <- ifelse(x[i]>k,muS[i],muF[i])   #change point process model

    p[i] ~ dnorm(mu[i],p.PC) #PhenoCam Data Model
    mn[i] ~ dnorm(mu[i],p.MN) # MODIS NDVI Data Model
    me[i] ~ dnorm(mu[i],p.ME) # MODIS EVI Data Model
    }
    }
    "
  } else{
    print("unknown season order")
}
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          n.chains = nchain)
  return(j.model)
  }
