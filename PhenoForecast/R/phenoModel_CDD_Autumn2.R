##' Creates a logistic phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @export
##' @import rjags
##' @import coda
phenoModel_CDD_Autumn2 <- function(data,nchain){
  ##Set priors
  data$s1 <- 0.5
  data$s2 <- 0.2
  data$fallLength.mu <- 450 ##Based on Richardson et al. (2007)
  data$fallLength.prec <- 1/(150**2) ##Based on Richardson et al. (2007)
  data$MOF.mu <- 400 ##CDD for middle of fall
  data$MOF.prec <- 1/(100**2)
  #data$sSlope.a <-
  #data$sSlope.b <-
  data$baseTemp <- 20


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
    for(i in 2:n){
      ##Calculate CDD
      Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])

      offset[i,yr] <- (baseTemp - Tair[i,yr])
      newCDD1[i,yr] <- ifelse(offset[i,yr]>0,offset[i,yr],0)

      newCDD2[i,yr] <- CDDs[(i-1),yr] + newCDD1[i,yr]
      CDDmu[i,yr] <- ifelse(Tair[i,yr]<baseTemp,newCDD2[i,yr],CDDs[(i-1),yr])
      CDDs[i,yr] ~ dnorm(CDDmu[i,yr],CDDprec)

      ##Calculate slopes for each time stamp
      m1l[i,yr] <- CDDs[i,yr]* CDDratio +b0_1
      m2l[i,yr] <- CDDs[i,yr]* -1* CDDratio +b0_2
      m1[i,yr] <- ifelse(CDDs[i,yr]<SOF,0,m1l)
      m2[i,yr] <- ifelse(CDDs[i,yr]<EOF,m2l,0)
      ml[i,yr] <- ifelse(CDDs[i,yr]<MOF,m1,m2)
      m[i,yr] <- min(ml,0)
      xl[i,yr] <- x[(i-1),yr] + m[i,yr]
      xl2[i,yr] ~ dnorm(xl[i,yr],p.proc)  ## process error
      x[i,yr] <- max(0, min(1,xl2[i,yr]) )
      }
  }


  for(i in 2:q){ ##Done for the current year forecast. Excluded from previous because n != q
    ##Calculate CDD
    Tair[i,N] ~ dnorm(TairMu[i,N],TairPrec[i,N])

    offset[i,N] <- (baseTemp - Tair[i,N])
    newCDD1[i,N] <- ifelse(offset[i,N]>0,offset[i,N],0)

    newCDD2[i,N] <- CDDs[(i-1),N] + newCDD1[i,N]
    CDDmu[i,N] <- ifelse(Tair[i,N]<baseTemp,newCDD2[i,N],CDDs[(i-1),N])
    CDDs[i,N] ~ dnorm(CDDmu[i,N],CDDprec)

    ##Calculate slopes for each time stamp
    m1l[i,N] <- CDDs[i,N]* CDDratio +b0_1
    m2l[i,N] <- CDDs[i,N]* -1* CDDratio +b0_2
    m1[i,N] <- ifelse(CDDs[i,N]<SOF,0,m1l)
    m2[i,N] <- ifelse(CDDs[i,N]<EOF,m2l,0)
    ml[i,N] <- ifelse(CDDs[i,N]<MOF,m1,m2)
    m[i,N] <- min(ml,0)
    xl[i,N] <- x[(i-1),N] + m[i,N]
    xl2[i,N] ~ dnorm(xl[i,N],p.proc)  ## process error
    x[i,N] <- max(0, min(1,xl2[i,N]) )

  }


  #### Priors
  for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dnorm(x_ic,tau_ic)
    CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
  }
  p.PC ~ dgamma(s1,s2)
  p.ME ~ dgamma(s1,s2)
  p.MN ~ dgamma(s2,s2)
  p.proc ~ dgamma(s1,s2)
  CDDprec ~ dgamma(s1,s2)
  fallLength ~ dnorm(fallLength.mu,fallLength.prec)
  MOF ~ dnorm(MOF.mu,MOF.prec)
  #sSlope ~ rbeta(sSlope.a,sSlope.b)
  sSlope ~ dunif(-1,0)

  ####Knowns based off of priors
  SOF <- MOF - fallLength/2
  EOF <- MOF + fallLength/2
  CDDratio <- (sSlope-0)/(MOF-SOF)
  b0_1 <- -1 * CDDratio * SOF
  b0_2 <- -1 * (-1*CDDratio) * EOF

  }"

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
