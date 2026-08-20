##' Creates a random walk phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, and data$n
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
randomWalkPhenoModel2020 <- function(data,nchain,priorCal=FALSE){
  ##Set priors
  if(typeof(priorCal)==typeof(FALSE)){
    data$s1.proc <- 1262.626
    data$s2.proc <- 50.50505
    data$x1.a <- 30 #Done to keep distribution close to 1 (over 75% of the data >0.95)
    data$x1.b <- 1
    data$s1.PC <- 1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
    data$s2.PC <- 50.50505 ##From mean <- 1/(0.2**2) and var = (mean-1/((0.4/1.96)**2))/2
  }

  ###JAGS model
  RandomWalk = "
    model{

    #### Data Models
    for(yr in 1:(N-1)){
      for(i in 1:n){
        p[i,yr] ~ dnorm(x[i,yr],p.PC)
      }
    }

    #### Process Model
    for(yr in 1:(N-1)){
      for(i in 2:n){
        xl[i,yr]~dnorm(x[i-1,yr],p.proc)
        x[i,yr] <- max(0, min(1,xl[i,yr]))
      }
    }

    #### Priors
    for(yr in 1:N){ ##Initial Conditions
      x[1,yr] ~ dbeta(x1.a,x1.b)
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)

    }
    "

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           n.chains = nchain)
  return(j.model)

  }
