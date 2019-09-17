##' Create a Bayes Model for a deciduous broadleaf site
##'
##' @param dataSource data source (GOES.NDVI, MODIS.NDVI, PC.GCC)
##' @param siteName Site Name
##' @param URL PhenoCam network URL
##' @param lat latitude of site in degrees
##' @param long longitude of site in degrees
##' @param startDay the day of year since 2017-01-01 to start the model
##' @param endDay the day of year since 2017-01-01 to end the model
##' @param niter the maximum number of iterations you want to give the model to converge within
##' @param seasonOrder
##' @import rjags
##' @import runjags
##' @export
createBayesModel.DB <- function(dataSource,siteName="",URL="",niter=100000,startDay,endDay,lat,long,TZ=5,PFT="DB",maxValue=FALSE,seasonOrder="FS") {
  nchain = 5
  inits <- list()
  if(dataSource=="PC.GCC"){
    data <- PC_data(siteName=siteName,URL=URL,startDate=as.Date("2017-07-01"),endDate=as.Date("2018-06-30"))
    #data$x[data$x<181] <- data$x[data$x<181]+365
    # inits.mu <- createInits(data=data,PFT=PFT)
    # for(i in 1:nchain){
    #   inits[[i]] <- list(TranS=rnorm(1,480,10),bS=rnorm(1,-0.10,0.015),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),k=rnorm(1,365,10))
    # }
    data$mean.c <- 0.3
    data$mean.d <- 0.35
    data$p.c <- 1/(0.15**2)
    data$p.d <- 1/(0.15**2)
    #print(range(data$x))
  }else if(dataSource == "MODIS.NDVI"){
    data = MODIS_data(siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay,metric="NDVI")
    print(length(data$x))
    data$obs.prec <- rep(1/((1.96/0.01)**2),length(data$x)) ##From Miura et al. (2000)
    print(data$obs.prec)
    data$mean.c <- 0.4
    data$p.c <- 1/(0.2**2)
    data$mean.d <- 0.6
    data$p.d <- 1/(0.2**2)
  }else if(dataSource == "MODIS.EVI"){
    data = MODIS_data(siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay,metric="EVI")

    data$mean.c <- 0.4
    data$obs.prec <- rep(1/((1.96/0.02)**2),length(data$x)) ##From Miura et al. (2000)
    print(data$obs.prec)
    data$p.c <- 1/(0.2**2)
    data$mean.d <- 0.6
    data$p.d <- 1/(0.2**2)
  }
  else if(dataSource=="GOES.NDVI"){
    if(maxValue){
      #print("Max values")
      data = GOES_data(siteName,startDay = startDay,endDay = endDay,lat=lat,long=long,TZ=TZ,maxValue=TRUE)
      #plot(data$x,data$y,pch=20)
    }
    else{
      data = GOES_data(siteName,startDay = startDay,endDay = endDay,lat=lat,long=long,TZ=TZ)
    }
    inits.mu <- createInits(data=data,PFT=PFT)
    for(i in 1:(nchain)){
      inits[[i]] <- list(TranS=rnorm(1,480,10),bS=rnorm(1,-0.09,0.05),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),k=rnorm(1,365,10))
    }
    data$mean.c <- 0.4
    data$p.c <- 1/(0.2**2)
    data$mean.d <- 0.6
    data$p.d <- 1/(0.2**2)
  }
  data$s1 <- 0.001
  data$s2 <- 0.00001
  data$p.Tran <- 1/(40**2)
  data$p.b <- 1/(0.05**2)
  data$mean.TranF <- 300
  data$mean.bF <- 0.10
  data$mean.TranS <- 475
  data$mean.bS <- -0.10
  data$mean.k <- 365
  data$p.k <- 1/(1**2)

  if(seasonOrder=="FS"){
    DB_model <- "
    model{
    ##priors
    TranS ~ dnorm(mean.TranS,p.Tran) ##S for spring
    bS ~ dnorm(mean.bS,p.b)
    TranF ~ dnorm(mean.TranF,p.Tran)  ##F for fall/autumn
    bF ~ dnorm(mean.bF,p.b)
    d ~ dnorm(mean.d,p.d)
    c ~ dnorm(mean.c,p.c)
    k ~ dnorm(mean.k,p.k)
    prec ~ dgamma(s1,s2)

    for(i in 1:n){
      muF[i] <- c/(1+exp(bF*(x[i]-TranF)))+d ##process model for fall
      muS[i] <- c/(1+exp(bS*(x[i]-TranS)))+d ##process model for Spring
      mu[i] <- ifelse(x[i]>k,muS[i],muF[i])   #change point process model

    y[i] ~ dnorm(mu[i],prec)
    yobs[i] ~ dnorm(y[i],obs.prec[i])
    }
    }
    "
  }
  else if (seasonOrder=="SF"){
    DB_model <- "
    model{
    ##priors
    TranS ~ dnorm(mean.TranS,p.Tran) ##S for spring
    bS ~ dnorm(mean.bS,p.b)
    TranF ~ dnorm(mean.TranF,p.Tran)  ##F for fall/autumn
    bF ~ dnorm(mean.bF,p.b)
    d ~ dnorm(mean.d,p.d)
    c ~ dnorm(mean.c,p.c)
    k ~ dnorm(mean.k,p.k)
    prec ~ dgamma(s1,s2)

    for(i in 1:n){
    muF[i] <- c/(1+exp(bF*(x[i]-TranF)))+d ##process model for fall
    muS[i] <- c/(1+exp(bS*(x[i]-TranS)))+d ##process model for Spring
    mu[i] <- ifelse(x[i]>k,muF[i],muS[i])   #change point process model

    y[i] ~ dnorm(mu[i],prec)
    yobs[i] ~ dnorm(y[i],obs.prec[i])
    }
    }
    "
  }
  else{
    print("unknown season order")
  }
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          n.chains = nchain)
  return(j.model)
}
