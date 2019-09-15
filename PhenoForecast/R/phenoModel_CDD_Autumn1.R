##' Creates a logistic phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
phenoModel_CDD_Autumn1 <- function(data,nchain){
  ##Set priors
  data$s1 <- 0.5
  data$s2 <- 0.2
  data$mu.b1 <- 0.2 #Based off of slope with points (sf=177,r=0) and (sf=250 and r = 1.5)
  data$prec.b1 <- 1/(0.05**2)
  data$mu.baseTemp <- 25 ##Kind of based off of Richardson et al. (2006)
  data$p.baseTemp <- 1/(2.5**2)
  data$mu.trans <- 600
  data$p.trans <- 1/(200**2)
  #data$mu.b0 <- -3.625 #Based off of slope with points (sf=177,r=0) and (sf=250 and r = 1.5)
  #data$prec.b0 <- 1/(0.5**2)

  ##JAGS code
  LogisticModel = "
  model{
  ### Data Models for complete years
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
      r[i,yr] <- ifelse(rl[i,yr]>r[(i-1),yr],r[(i-1),yr],rl[i,yr])
    }
    for(i in 2:n){
      Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
      SfNew[i,yr] <- Sf[(i-1),yr] + (baseTemp - Tair[i,yr])
      Sfmu[i,yr] <- ifelse(Tair[i,yr]<baseTemp,SfNew[i,yr],Sf[(i-1),yr])
      Sf[i,yr] ~ dnorm(Sfmu[i,yr],Sfprec)

      rl[i,yr] <- b1 * Sf[i,yr] + b0
      color[i,yr] <- x[(i-1),yr] + r[i,yr] * x[(i-1),yr] * (1-x[(i-1),yr])  ## latent process

      #Sf[i,yr] ~ dnorm(Sfmu[i,yr],Sfprec[i,yr])
      xl[i,yr] ~ dnorm(color[i,yr],p.proc)  ## process error
      x[i,yr] <- max(0, min(1,xl[i,yr]) ) ## trunate normal process error
    }
  }
  r[2,N] <- rl[2,N]
  for(i in 3:q){
    r[i,N] <- ifelse(rl[i,N]>r[(i-1),N],r[(i-1),N],rl[i,N])
  }

  for(i in 2:q){ ##Done for the current year forecast. Excluded from previous because n != q
    Tair[i,N] <- dnorm(TairMu[i,N],TairPrec[i,N])
    SfNew[i,N] <- Sf[(i-1),N] + (baseTemp - Tair[i,N])
    Sfmu[i,N] <- ifelse(Tair[i,N]<baseTemp,SfNew[i,N],Sf[(i-1),N])
    Sf[i,N] ~ dnorm(Sfmu[i,N],Sfprec)

    rl[i,N] <- b1 * Sf[i,N] + b0
    color[i,N] <- x[(i-1),N] + r[i,N] * x[(i-1),N] * (1-x[(i-1),N])  ## latent process

    xl[i,N] ~ dnorm(color[i,N],p.proc)  ## process error
    x[i,N] <- max(0, min(1,xl[i,N]) ) ## trunate normal process error
  }

  b0 <-  b1 * -1 * trans

  #### Priors
  for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dnorm(x_ic,tau_ic)
    r[1,yr] ~ dnorm(-0.4,100)
    color[1,yr] ~ dnorm(x_ic,tau_ic)
    colorT[1,yr] ~ dnorm(x_ic,tau_ic)
    Sf[1,N] <- ifelse(Tair[1,yr]<baseTemp,Tair[1,yr],0)
  }
  p.PC ~ dgamma(s1,s2)
  p.ME ~ dgamma(s1,s2)
  p.MN ~ dgamma(s2,s2)
  p.proc ~ dgamma(s1,s2)
  #trans ~ dnorm(110,0.000625)
  trans ~ dnorm(mu.trans,p.trans)
  Sfprec ~ dgamma(s1,s2)

  b1 ~ dnorm(mu.b1,prec.b1)
  baseTemp ~ dnorm(mu.baseTemp,p.baseTemp)
  }"

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
