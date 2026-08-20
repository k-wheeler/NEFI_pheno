##' Creates a random walk phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, and data$n
##' @param nchain The desired number of chains in the MCMC
##' @param index The desired index (all, GCC, NDVI, or EVI)
##' @export
##' @import rjags
##' @import coda
randomWalkPhenoModel <- function(data,nchain,index="all"){
  ##Set priors

  data$s1.proc <- 1262.626
  data$s2.proc <- 50.50505
  data$x1.a <- 1 #Done to keep distribution close to 0 (over 75% of the data < 0.05)
  data$x1.b <- 30
  if(index=="GCC"){
    data$mn <- NULL
    data$me <- NULL
    data$s1.PC <- 1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
    data$s2.PC <- 50.50505 ##From mean <- 1/(0.2**2) and var = (mean-1/((0.4/1.96)**2))/2

    ###JAGS model
    RandomWalk = "
    model{

    #### Data Models
    for(i in 1:n){
    p[i] ~ dnorm(x[i],p.PC) #PhenoCam Data Model
    #mn[i] ~ dnorm(x[i],p.MN) # MODIS NDVI Data Model
    #me[i] ~ dnorm(x[i],p.ME) # MODIS EVI Data Model
    }

    #### Process Model
    for(i in 2:n){
    xl[i]~dnorm(x[i-1],p.proc)
    x[i] <- max(0, min(1,xl[i]))
    }

    #### Priors
    x[1] ~ dbeta(x1.a,x1.b)
    p.PC ~ dgamma(s1.PC,s2.PC)
    #p.MN ~ dgamma(s1.MN,s2.MN)
    #p.ME ~ dgamma(s1.ME,s2.ME)
    p.proc ~ dgamma(s1.proc,s2.proc)

    }
    "

  }else{
    data$s1.PC <- 1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
    data$s2.PC <- 50.50505 ##From mean <- 1/(0.2**2) and var = (mean-1/((0.4/1.96)**2))/2
    data$s1.MN <- 76880.05 ##From mean = 1/((0.01/1.96)**2) and var = (mean - 1/((0.4/1.96)**2))/2
    data$s2.MN <- 2.001251
    data$s1.ME <- 19256.14 ##From mean = 1/((0.02/1.96)**2) and var = (mean - 1/((0.4/1.96)**2))/2
    data$s2.ME <- 2.005013



  ###JAGS model
  RandomWalk = "
  model{

  #### Data Models
  for(i in 1:n){
    p[i] ~ dnorm(x[i],p.PC) #PhenoCam Data Model
    mn[i] ~ dnorm(x[i],p.MN) # MODIS NDVI Data Model
    me[i] ~ dnorm(x[i],p.ME) # MODIS EVI Data Model
  }

  #### Process Model
  for(i in 2:n){
    xl[i]~dnorm(x[i-1],p.proc)
    x[i] <- max(0, min(1,xl[i]))
  }

  #### Priors
  x[1] ~ dbeta(x1.a,x1.b)
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.MN ~ dgamma(s1.MN,s2.MN)
  p.ME ~ dgamma(s1.ME,s2.ME)
  p.proc ~ dgamma(s1.proc,s2.proc)

  }
  "
  }

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
