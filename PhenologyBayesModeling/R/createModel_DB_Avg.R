##' Create a Bayes Model for a deciduous broadleaf site using the window average GOES values
##'
##' @param dataSource data source (GOES.NDVI, MODIS.NDVI, PC.GCC)
##' @param siteName Site Name
##' @param URL PhenoCam network URL
##' @param lat latitude of site in degrees
##' @param long longitude of site in degrees
##' @param startDay the day of year since 2017-01-01 to start the model
##' @param endDay the day of year since 2017-01-01 to end the model
##' @param niter the maximum number of iterations you want to give the model to converge within
##' @import rjags
##' @import runjags
##' @export
createBayesModel.DB_Avg <- function(dataSource="GOES.NDVI",siteName="",URL="",niter=100000,startDay,endDay,lat,long,TZ=5) {
  nchain = 5
  inits <- list()
  # if(dataSource=="PC.GCC"){
  #   data <- PC_data(siteName=siteName,URL=URL,startDay=startDay,endDay=endDay)
  #   inits.mu <- createInits(data=data,PFT=PFT)
  #   for(i in 1:nchain){
  #     inits[[i]] <- list(TranS=rnorm(1,480,10),bS=rnorm(1,-0.10,0.015),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),k=rnorm(1,365,10))
  #   }
  #   data$mean.c <- 0.1
  #   data$mean.d <- 0.35
  #   data$p.c <- 1/(0.1**2)
  #   data$p.d <- 1/(0.1**2)
  # }
  # else if(dataSource == "MODIS.NDVI"){
  #   data = MODIS_data(siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay,metric="NDVI")
  #   inits.mu <- createInits(data=data,PFT=PFT)
  #   for(i in 1:(nchain)){
  #     inits[[i]] <- list(TranS=rnorm(1,480,10),bS=rnorm(1,-0.09,0.05),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),k=rnorm(1,365,10))
  #   }
  #   data$mean.c <- 0.4
  #   data$p.c <- 1/(0.2**2)
  #   data$mean.d <- 0.6
  #   data$p.d <- 1/(0.2**2)
  # }
  # else if(dataSource == "MODIS.EVI"){
  #   data = MODIS_data(siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay,metric="EVI")
  #   inits.mu <- createInits(data=data,PFT=PFT)
  #   for(i in 1:(nchain)){
  #     inits[[i]] <- list(TranS=rnorm(1,480,10),bS=rnorm(1,-0.09,0.05),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),k=rnorm(1,365,10))
  #   }
  #   data$mean.c <- 0.4
  #   data$p.c <- 1/(0.2**2)
  #   data$mean.d <- 0.6
  #   data$p.d <- 1/(0.2**2)
  # }
  if(dataSource=="GOES.NDVI"){
    data = GOES_data(siteName,startDay = startDay,endDay = endDay,lat=lat,long=long,TZ=TZ,window=TRUE)
    #inits.mu <- createInits(data=data,PFT="DB")
    for(i in 1:(nchain)){
      #inits[[i]] <- list(TranS=rnorm(1,480,10),bS=rnorm(1,-0.09,0.05),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),k=rnorm(1,365,10))
      inits[[i]] <- list(TranS=rnorm(1,480,10),bS=rnorm(1,-0.09,0.05),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,0.4,0.02),d=rnorm(1,0.3,0.02),k=rnorm(1,365,10))

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
  data$p.k <- 1/(30**2)
  #print(data$n)

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
  prec ~ dgamma(s1,s2) #process error

  for(i in 1:n){
  muF[i] <- c/(1+exp(bF*(x[i]-TranF)))+d ##process model for fall
  muS[i] <- c/(1+exp(bS*(x[i]-TranS)))+d ##process model for Spring
  mu[i] <- ifelse(x[i]>k,muS[i],muF[i])   #change point process model

  y[i] ~ dnorm(mu[i],prec)
  w[i] <- obs.prec[i] * size[i]
  yobs[i] ~ dnorm(y[i],w[i])
  }
  }
  "
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          inits=inits,
                          n.chains = nchain)

  return(j.model)
}

#y[i]  ~ dnorm(mu[i],prec[i])		## data model for logistic (need to divide betwen obs and process error or could add them together because they are precisions)
#having them separate vs together might have different timing
#y[i] ~ dnorm(mu[i],prec+prec.obs[i])
#want to keep it as Gibbs and not MH (look at trace plots to see which converges quicker)(can look at acceptance rate of prec to see if it switches to another sampler)
#might be slower because it has a extra dummy variable that you don't want to spit out
