##' Creates a logistic phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
logisticCovPhenoModel <- function(data,nchain){
  ##Set priors
  data$s1 <- 0.5
  data$s2 <- 0.2
  data$min.b0 <- -1000
  data$max.b0 <- 0 ##Intercept as always negative
  data$min.b1 <- 0 ##Slope is always positive
  data$max.b1 <- 1000

  ##JAGS code
  LogisticModel = "
  model{
  #### Data Models for complete years
  for(yr in 1:(N-1)){
    for(i in 1:n){
      p[i,yr] ~ dnorm(x[i,yr],p.PC)
      mn[i,yr] ~ dnorm(x[i,yr],p.MN)
      me[i,yr] ~ dnorm(x[i,yr],p.ME)
    }
  }
  ##Data Model for last year
  for(i in 1:q){
      p[i,N] ~ dnorm(x[i,N],p.PC)
      mn[i,N] ~ dnorm(x[i,N],p.MN)
      me[i,N] ~ dnorm(x[i,N],p.ME)
  }

  #### Process Model
  for(yr in 1:(N-1)){
    for(i in 2:n){
      r[i,yr] <- b0 + b1 * Sf[i,yr]
      color[i,yr] <- x[(i-1),yr] + r[i,yr] * x[(i-1),yr] * (1-x[(i-1),yr])  ## latent process
      Sf[i,yr] ~ dnorm(Sfmu[i,yr],Sfprec[i,yr])
      xl[i,yr] ~ dnorm(color[i,yr],p.proc)  ## process error
      x[i,yr] <- max(0, min(1,xl[i,yr]) ) ## trunate normal process error
    }
  }
  for(i in 2:q){ ##Done for the current year forecast. Excluded from previous because n != q
      r[i,N] <- b0 + b1 * Sf[i,N]
      color[i,N] <- x[(i-1),N] + r[i,N] * x[(i-1),N] * (1-x[(i-1),N])  ## latent process
      Sf[i,N] ~ dnorm(Sfmu[i,N],Sfprec[i,N])
      xl[i,N] ~ dnorm(color[i,N],p.proc)  ## process error
      x[i,N] <- max(0, min(1,xl[i,N]) ) ## trunate normal process error
  }

  #### Priors
  for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dnorm(x_ic,tau_ic)
    color[1,yr] ~ dnorm(x_ic,tau_ic)
  }
  p.PC ~ dgamma(s1,s2)
  p.ME ~ dgamma(s1,s2)
  p.MN ~ dgamma(s2,s2)
  p.proc ~ dgamma(s1,s2)
  b0 ~ dunif(min.b0,max.b0) ##Need to change to normal
  b1 ~ dunif(min.b1,max.b1)
  }"

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
