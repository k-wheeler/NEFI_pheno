##' Creates a logistic phenology forecast model based on PhenoCam
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @param baseTemp The desired base temperature
##' @export
##' @import rjags
##' @import coda
phenoModel_CDD_Autumn_GCC <- function(data,nchain,MODIS_index="NDVI",baseTemp=NA){
  inits <- list()
  for(i in 1:nchain){
    inits[[i]] <- list(p.PC=rnorm(1,30,0.2),p.proc=rnorm(1,26,0.2),
                       MOF=rnorm(1,500,50),fallLength=rnorm(1,500,100),
                       sSlope=rnorm(1,0.05,0.002))
  }
  ##Set priors
  data$s1.PC <- 8.4##1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
  data$s2.PC <- 0.02##50.50505 ##From mean <- 1/(0.2**2) and var = (mean-1/((0.4/1.96)**2))/2
  #data$s1.MN <- 76880.05 ##From mean = 1/((0.01/1.96)**2) and var = (mean - 1/((0.4/1.96)**2))/2
  #data$s2.MN <- 2.001251
  #data$s1.ME <- 19256.14 ##From mean = 1/((0.02/1.96)**2) and var = (mean - 1/((0.4/1.96)**2))/2
  #data$s2.ME <- 2.005013
  data$s1.proc <- 8.4##1262.626
  data$s2.proc <- 0.02##50.50505
  data$x1.a <- 30 #Done to keep distribution close to 0 (over 75% of the data < 0.05)
  data$x1.b <- 1
  data$alp.sSlope <- 2.5 ##Based on investigating O'Keefe Harvard phenology data
  data$bet.sSlope <- 23 ##Based on investigating O'Keefe Harvard phenology data

  data$fallLength.mu <- 450 ##Based on Richardson et al. (2006)
  data$fallLength.prec <- 1/(150**2) ##Based on Richardson et al. (2006)
  data$MOF.mu <- 400 ##CDD for middle of fall
  data$MOF.prec <- 1/(100**2)

  ##JAGS code
  if(is.na(baseTemp)){
    data$baseTemp.mu <- 20
    data$baseTemp.prec <- 1/(10**2)
    LogisticModel = "
    model{
    ### Data Models for complete years
    for(yr in 1:(N-1)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    ##Data Model for current year
    for(i in 1:q){
    p[i,N] ~ dnorm(x[i,N],p.PC)
    }

    #### Process Model
    for(yr in 1:(N-1)){
    for(i in 2:n){
    ##Calculate CDD
    Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])

    offset[i,yr] <- (baseTemp - Tair[i,yr])
    newCDD1[i,yr] <- ifelse(offset[i,yr]>0,offset[i,yr],0)

    newCDD2[i,yr] <- CDDs[(i-1),yr] + newCDD1[i,yr]
    CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,newCDD2[i,yr],CDDs[(i-1),yr])
    #CDDs[i,yr] ~ dnorm(CDDmu[i,yr],CDDprec)

    ##Calculate slopes for each time stamp
    m1l[i,yr] <- CDDs[i,yr] * CDDratio + b0_1
    m2l[i,yr] <- CDDs[i,yr] * -1* CDDratio + b0_2
    m1[i,yr] <- ifelse(CDDs[i,yr]<SOF,0,m1l[i,yr])
    m2[i,yr] <- ifelse(CDDs[i,yr]<EOF,m2l[i,yr],0)
    ml[i,yr] <- ifelse(CDDs[i,yr]<MOF,m1[i,yr],m2[i,yr])
    m[i,yr] <- min(ml[i,yr],0)
    xl[i,yr] <- x[(i-1),yr] + m[i,yr]
    xl2[i,yr] ~ dnorm(xl[i,yr],p.proc)  ## process error
    x[i,yr] <- max(0.01, min(1,xl2[i,yr]) )
    }
    }


    for(i in 2:q){ ##Done for the current year forecast. Excluded from previous because n != q
    ##Calculate CDD
    Tair[i,N] ~ dnorm(TairMu[i,N],TairPrec[i,N])

    offset[i,N] <- (baseTemp - Tair[i,N])
    newCDD1[i,N] <- ifelse(offset[i,N]>0,offset[i,N],0)

    newCDD2[i,N] <- CDDs[(i-1),N] + newCDD1[i,N]
    CDDs[i,N] <- ifelse(Tair[i,N]<baseTemp,newCDD2[i,N],CDDs[(i-1),N])
    #CDDs[i,N] ~ dnorm(CDDmu[i,N],CDDprec)

    ##Calculate slopes for each time stamp
    m1l[i,N] <- CDDs[i,N] * CDDratio + b0_1
    m2l[i,N] <- CDDs[i,N] * -1 * CDDratio + b0_2
    m1[i,N] <- ifelse(CDDs[i,N]<SOF,0,m1l[i,N])
    m2[i,N] <- ifelse(CDDs[i,N]<EOF,m2l[i,N],0)
    ml[i,N] <- ifelse(CDDs[i,N]<MOF,m1[i,N],m2[i,N])
    m[i,N] <- min(ml[i,N],0)
    xl[i,N] <- x[(i-1),N] + m[i,N]
    xl2[i,N] ~ dnorm(xl[i,N],p.proc)  ## process error
    x[i,N] <- max(0, min(1,xl2[i,N]) )

    }


    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a,x1.b)
    CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    #p.MN ~ dgamma(s1.MN,s2.MN)
    #p.ME ~ dgamma(s1.ME,s2.ME)
    p.proc ~ dgamma(s1.proc,s2.proc)

    fallLength ~ dnorm(fallLength.mu,fallLength.prec)
    MOF ~ dnorm(MOF.mu,MOF.prec)
    sSlope ~ dunif(-1,0)
    baseTemp ~ dnorm(baseTemp.mu,baseTemp.prec)

    ####Knowns based off of priors
    SOF <- MOF - fallLength/2
    EOF <- MOF + fallLength/2
    CDDratio <- (sSlope-0)/(MOF-SOF)
    b0_1 <- -1 * CDDratio * SOF
    b0_2 <- -1 * (-1*CDDratio) * EOF

    }"
  }else{
    data$baseTemp <- baseTemp
    LogisticModel = "
    model{
    ### Data Models for complete years
    for(yr in 1:(N-1)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    ##Data Model for current year
    for(i in 1:q){
    p[i,N] ~ dnorm(x[i,N],p.PC)
    }

    #### Process Model
    for(yr in 1:(N-1)){
    for(i in 2:n){
    ##Calculate CDD
    Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])

    offset[i,yr] <- (baseTemp - Tair[i,yr])
    newCDD1[i,yr] <- ifelse(offset[i,yr]>0,offset[i,yr],0)

    newCDD2[i,yr] <- CDDs[(i-1),yr] + newCDD1[i,yr]
    CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,newCDD2[i,yr],CDDs[(i-1),yr])
    #CDDs[i,yr] ~ dnorm(CDDmu[i,yr],CDDprec)

    ##Calculate slopes for each time stamp
    m1l[i,yr] <- CDDs[i,yr] * CDDratio + b0_1
    m2l[i,yr] <- CDDs[i,yr] * -1* CDDratio + b0_2
    m1[i,yr] <- ifelse(CDDs[i,yr]<SOF,0,m1l[i,yr])
    m2[i,yr] <- ifelse(CDDs[i,yr]<EOF,m2l[i,yr],0)
    ml[i,yr] <- ifelse(CDDs[i,yr]<MOF,m1[i,yr],m2[i,yr])
    xl[i,yr] <- x[(i-1),yr] + min(ml[i,yr],0)
    #xl[i,yr] <- x[(i-1),yr] + m[i,yr]
    #is_censored[i,yr] ~ dinterval(p[i,yr], c(0,x[(i-1),yr]))
    xl2[i,yr] ~ dnorm(xl[i,yr],p.proc)  ## process error
    #xl3[i,yr] <- min(xl[(i-1),yr],xl2[i,yr])
    x[i,yr] <- max(0, min(1,xl2[i,yr]) )
    }
    }


    for(i in 2:q){ ##Done for the current year forecast. Excluded from previous because n != q
    ##Calculate CDD
    Tair[i,N] ~ dnorm(TairMu[i,N],TairPrec[i,N])

    offset[i,N] <- (baseTemp - Tair[i,N])
    newCDD1[i,N] <- ifelse(offset[i,N]>0,offset[i,N],0)

    newCDD2[i,N] <- CDDs[(i-1),N] + newCDD1[i,N]
    CDDs[i,N] <- ifelse(Tair[i,N]<baseTemp,newCDD2[i,N],CDDs[(i-1),N])
    #CDDs[i,N] ~ dnorm(CDDmu[i,N],CDDprec)

    ##Calculate slopes for each time stamp
    m1l[i,N] <- CDDs[i,N] * CDDratio + b0_1
    m2l[i,N] <- CDDs[i,N] * -1 * CDDratio + b0_2
    m1[i,N] <- ifelse(CDDs[i,N]<SOF,0,m1l[i,N])
    m2[i,N] <- ifelse(CDDs[i,N]<EOF,m2l[i,N],0)
    ml[i,N] <- ifelse(CDDs[i,N]<MOF,m1[i,N],m2[i,N])
    xl[i,N] <- x[(i-1),N] + min(ml[i,N],0)
    xl2[i,N] ~ dnorm(xl[i,N],p.proc)  ## process error
    #is_censored[i,N] ~ dinterval(p[i,N], c(0,x[(i-1),N]))
    #xl3[i,N] <- min(xl[(i-1),N],xl2[i,N])
    x[i,N] <- max(0, min(1,xl2[i,N]) )

    }


    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a,x1.b)
    xl[1,yr] ~ dbeta(x1.a,x1.b)
    xl3[1,yr] ~ dbeta(x1.a,x1.b)
    CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    #p.MN ~ dgamma(s1.MN,s2.MN)
    #p.ME ~ dgamma(s1.ME,s2.ME)
    p.proc ~ dgamma(s1.proc,s2.proc)

    fallLength ~ dnorm(fallLength.mu,fallLength.prec)
    MOF ~ dnorm(MOF.mu,MOF.prec)
    sSlope ~ dbeta(alp.sSlope,bet.sSlope)
    sSlopeNeg <- sSlope * -1
    #sSlope ~ dunif(-1,0) ##Steepest slope
    #baseTemp ~ dnorm(baseTemp.mu,baseTemp.prec)

    ####Knowns based off of priors
    SOF <- MOF - fallLength/2
    EOF <- MOF + fallLength/2
    CDDratio <- (sSlopeNeg-0)/(MOF-SOF)
    b0_1 <- -1 * CDDratio * SOF
    b0_2 <- -1 * (-1*CDDratio) * EOF

  }"
  }

  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           inits = inits,
                           n.chains = nchain)
  return(j.model)

}
