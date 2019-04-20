##' Creates a logistic phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
logisticCovPhenoModel2 <- function(data,nchain){
  ##Set priors
  data$s1 <- 0.5
  data$s2 <- 0.2
  data$mu.b1 <- 0.2 #Based off of slope with points (sf=177,r=0) and (sf=250 and r = 1.5)
  data$prec.b1 <- 1/(0.05**2)
  #data$mu.b0 <- -3.625 #Based off of slope with points (sf=177,r=0) and (sf=250 and r = 1.5)
  #data$prec.b0 <- 1/(0.5**2)

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
  ##Data Model for current year
  for(i in 1:q){
      p[i,N] ~ dnorm(x[i,N],p.PC)
      mn[i,N] ~ dnorm(x[i,N],p.MN)
      me[i,N] ~ dnorm(x[i,N],p.ME)
  }

  #### Process Model
  for(yr in 1:(N-1)){
    r[2,yr] <- rl[2,yr]
    for(i in 3:n){
      r[i,yr] <- ifelse(rl[i,yr]<r[(i-1),yr],rl[(i-1),yr],rl[i,yr])
    }
    for(i in 2:n){
      rl[i,yr] <- b1 * Sf[i,yr] + b0
      color[i,yr] <- x[(i-1),yr] + r[i,yr] * x[(i-1),yr] * (1-x[(i-1),yr])  ## latent process
      Sf[i,yr] ~ dnorm(Sfmu[i,yr],Sfprec[i,yr])
      xl[i,yr] ~ dnorm(color[i,yr],p.proc)  ## process error
      x[i,yr] <- max(0, min(1,xl[i,yr]) ) ## trunate normal process error
    }
  }
  r[2,N] <- rl[2,N]
  for(i in 3:q){
    r[i,N] <- ifelse(rl[i,N]<r[(i-1),N],rl[(i-1),N],rl[i,N])
  }
  for(i in 2:q){ ##Done for the current year forecast. Excluded from previous because n != q
      rl[i,N] <- b1 * Sf[i,N] + b0
      color[i,N] <- x[(i-1),N] + r[i,N] * x[(i-1),N] * (1-x[(i-1),N])  ## latent process
      Sf[i,N] ~ dnorm(Sfmu[i,N],Sfprec[i,N])
      xl[i,N] ~ dnorm(color[i,N],p.proc)  ## process error
      x[i,N] <- max(0, min(1,xl[i,N]) ) ## trunate normal process error
  }

  b0 <-  b1 * -1 * trans

  #### Priors
  for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dnorm(x_ic,tau_ic)
    r[1,yr] ~ dnorm(-0.4,100)
    color[1,yr] ~ dnorm(x_ic,tau_ic)
  }
  p.PC ~ dgamma(s1,s2)
  p.ME ~ dgamma(s1,s2)
  p.MN ~ dgamma(s2,s2)
  p.proc ~ dgamma(s1,s2)
  trans ~ dnorm(300,0.005)
  #b0 ~ dunif(min.b0,max.b0) ##Need to change to normal
  #b1 ~ dunif(min.b1,max.b1)
  b1 ~ dnorm(mu.b1,prec.b1)
 # b0 ~ dnorm(mu.b0,prec.b0)
  }"

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
