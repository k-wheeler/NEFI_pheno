##' Creates a random walk phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, and data$n
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
randomWalkPhenoModel <- function(data,nchain){
  ##Set priors (note: this could be condensed for now since all of the priors are the same)
  data$s1.PC <- 0.5
  data$s2.PC <- 0.2
  data$s1.MN <- 0.5
  data$s2.MN <- 0.2
  data$s1.ME <- 0.5
  data$s2.ME <- 0.2
  data$s1.proc <- 0.5
  data$s2.proc <- 0.2

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
  x[i]~dnorm(x[i-1],p.proc)
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.MN ~ dgamma(s1.MN,s2.MN)
  p.ME ~ dgamma(s1.ME,s2.ME)
  p.proc ~ dgamma(s1.proc,s2.proc)
  }
  "

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
