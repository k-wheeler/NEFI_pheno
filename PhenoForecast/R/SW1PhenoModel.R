##' Creates the first spring warming model from Melaas 2016 GCB
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
SW1PhenoModel <- function(data,nchain){

  ##Set priors (note: this is condensed for now since all of these priors are the same)
  data$s1 <- 0.5
  data$s2 <- 0.2

  print(dim(data$mn))
  print(dim(data$me))
  print(dim(data$p))
  print(data$N)
  print(data$n)

  ###JAGS model
  SW1= "
  model{

  #### Data Models
  for(yr in 1:N){
    for(i in 1:n){
      mn[i,yr] ~ dnorm(x[i,yr],p.MN) ##MODIS NDVI data model
      p[i,yr] ~ dnorm(x[i,yr],p.PC) ##PhenoCam data model
      me[i,yr] ~ dnorm(x[i,yr],p.ME) ##MODIS EVI data model
    }
  }

  #### Process Model
  # for(yr in 1:N){
  #   for(i in 2:n){
  #     mu.include[i,yr] <- (28.4)/(1+exp(3.4-0.185*T.air[i,yr]))
  #     mu.exclude[i,yr] <- 0
  #     R.f[i,yr] <- ifelse(T.air[i,yr]>0,mu.include[i,yr],mu.exclude[i,yr])
  #
  #     S.f[i,yr] <- S.f[(i-1),yr]+R.f[i,yr] ##Summation
  #
  #     #x[i,yr] ~ dnorm(x[(i-1),yr],p.proc)
  #     #F.stars[i,yr] <- ifelse(x[i,yr]==0.20,1,0) ##S.f[i,yr]
  #   }
  #   #F.star.yr[yr] <- max(F.stars[,yr])
  #   #sigma.yr[yr] <- F.star.yr[yr] - F.star
  # }


  #### Priors
  x[1,1] ~ dnorm(x_ic,tau_ic)
  x[1,2] ~ dnorm(x_ic,tau_ic)
  x[1,3] ~ dnorm(x_ic,tau_ic)
  x[1,4] ~ dnorm(x_ic,tau_ic)
  x[1,5] ~ dnorm(x_ic,tau_ic)
  x[1,6] ~ dnorm(x_ic,tau_ic)
  x[1,7] ~ dnorm(x_ic,tau_ic)
  x[1,8] ~ dnorm(x_ic,tau_ic)
  p.PC ~ dgamma(s1,s2)
  p.MN ~ dgamma(s1,s2)
  p.ME ~ dgamma(s1,s2)
  #p.proc ~ dgamma(s1.proc,s2.proc)
  p.F ~ dgamma(s1,s2)
  #F.star ~ unif(0,100000)
  }
  "

  ###Create the JAGS model

  j.model   <- jags.model (file = textConnection(SW1),
                           data = data,
                           n.chains = nchain)
  return(j.model)




}
