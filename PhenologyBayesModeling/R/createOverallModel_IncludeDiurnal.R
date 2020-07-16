##' Create a Bayes Model for a deciduous broadleaf site based on the diurnal-fit NDVI data
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
createBayesModel.DB_Overall <- function(data,niter=100000,seasonOrder="FS") {
  nchain = 5
  inits <- list()
  inits.mu <- createInits(data=data,PFT="DB")
  for(i in 1:(nchain)){
    inits[[i]] <- list(TranS=rnorm(1,115,10),bS=rnorm(1,-0.09,0.05),TranF=rnorm(1,280,10),bF=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),k=rnorm(1,365,10))
  }
  data$mean.c <- 0.4
  data$p.c <- 1/(0.2**2)
  data$mean.d <- 0.6
  data$p.d <- 1/(0.2**2)

  data$s1 <- 0.001
  data$s2 <- 0.00001
  data$p.Tran <- 1/(40**2)
  data$p.b <- 1/(0.05**2)
  data$mean.TranF <- 300
  data$mean.bF <- 0.10

  data$mean.bS <- -0.10

  #print(data$n)
  #print(data$x)
  #print(data$y)
  #print(data$obs.prec)

  if(seasonOrder=="FS"){
    data$mean.TranS <- 475
    data$mean.k <- 365
    data$p.k <- 1/(30**2)

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
  yobs[i] ~ dnorm(y[i],obs.prec[i])
  }
  }
  "
  }else if(seasonOrder=="SF"){
    data$mean.TranS <- 475-365
    data$k <- 182
    DB_model <- "
    model{
    ##priors
    TranS ~ dnorm(mean.TranS,p.Tran) ##S for spring
    bS ~ dnorm(mean.bS,p.b)
    TranF ~ dnorm(mean.TranF,p.Tran)  ##F for fall/autumn
    bF ~ dnorm(mean.bF,p.b)
    d ~ dnorm(mean.d,p.d)
    c ~ dnorm(mean.c,p.c)
    #k ~ dnorm(mean.k,p.k)
    prec ~ dgamma(s1,s2) #process error

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
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          inits=inits,
                          n.chains = nchain)

  return(j.model)
}
