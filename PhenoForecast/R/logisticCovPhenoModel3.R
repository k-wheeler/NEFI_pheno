##' Creates a logistic phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
logisticCovPhenoModel3 <- function(data,nchain){
  ##Set priors
  data$s1 <- 0.5
  data$s2 <- 0.2
  data$mu.b1 <- 0.2 #Based off of slope with points (sf=177,r=0) and (sf=250 and r = 1.5)
  data$prec.b1 <- 1/(0.05**2)

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
    for(i in 2:n){
      rl[i,yr] <- b1 * Sf[i,yr] + b0
      r[i,yr] <- ifelse(rl[i,yr]<0,0,rl[i,yr])
      color[i,yr] <- x[(i-1),yr] + r[i,yr] * x[(i-1),yr] * (1-x[(i-1),yr])  ## latent process
      Sf[i,yr] ~ dnorm(Sfmu[i,yr],Sfprec[i,yr])
      xl[i,yr] ~ dnorm(color[i,yr],p.proc)  ## process error
      x[i,yr] <- max(0, min(1,xl[i,yr]) ) ## trunate normal process error
    }
  }
  for(i in 2:q){ ##Done for the current year forecast. Excluded from previous because n != q
    rl[i,N] <- b1 * Sf[i,N] + b0
    r[i,N] <- ifelse(rl[i,N]<0,0,rl[i,N])
    color[i,N] <- x[(i-1),N] + r[i,N] * x[(i-1),N] * (1-x[(i-1),N])  ## latent process

    Sf[i,N] ~ dnorm(Sfmu[i,N],Sfprec[i,N])
    xl[i,N] ~ dnorm(color[i,N],p.proc)  ## process error
    x[i,N] <- max(0, min(1,xl[i,N]) ) ## trunate normal process error
  }

  b0 <-  b1 *-1 * trans

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
  #trans ~ dnorm(110,0.000625)I(0,182)
  trans ~ dnorm(300,0.005)
  b1 ~ dnorm(mu.b1,prec.b1)
  # b0 ~ dnorm(mu.b0,prec.b0)
  }"

  # inits <- list()
  # for(n in 1:nchain){
  #   inits[[n]] <- list(trans=rnorm(1,110,10),b1=rnorm(1,0.2,0.01))
  # }
  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}

