##' Creates a logistic phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @param season Season: spring or fall
##' @param index The desired index (all, GCC, NDVI, or EVI)
##' @export
##' @import rjags
##' @import coda
logisticPhenoModel <- function(data,nchain,season,index){
  ##Set priors

  if(index=="GCC"){
    data$mn <- NULL
    data$me <- NULL
    data$s1.PC <- 1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
    data$s2.PC <- 50.50505 ##From mean <- 1/(0.2**2) and var = (mean-1/((0.4/1.96)**2))/2
  }else{
    data$s1.PC <- 1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
    data$s2.PC <- 50.50505 ##From mean <- 1/(0.2**2) and var = (mean-1/((0.4/1.96)**2))/2
    data$s1.MN <- 76880.05 ##From mean = 1/((0.01/1.96)**2) and var = (mean - 1/((0.4/1.96)**2))/2
    data$s2.MN <- 2.001251
    data$s1.ME <- 19256.14 ##From mean = 1/((0.02/1.96)**2) and var = (mean - 1/((0.4/1.96)**2))/2
    data$s2.ME <- 2.005013
  }

  data$s1.proc <- 1262.626
  data$s2.proc <- 50.50505
  if(season=="spring"){
    data$x1.a <- 1 #Done to keep distribution close to 0 (over 75% of the data < 0.05)
    data$x1.b <- 30
  }else{
    data$x1.a <- 30 #Done to keep distribution close to 1 (over 75% of the data > 0.95)
    data$x1.b <- 1
  }


  ##JAGS code
  if(season=="spring"){
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
    #### Color is the expected new phenology stage given the previous stage and logistic
    #### subtraction instead of addition in the discrete logistic eqn makes r negative (so logistic goes down).
    for(yr in 1:(N-1)){
      for(i in 2:n){
        color[i,yr] <- x[(i-1),yr] + r * x[(i-1),yr] * (1-x[(i-1),yr])  ## latent process
        xl[i,yr] ~ dnorm(color[i,yr],p.proc)  ## process error
        x[i,yr] <- max(0, min(1,xl[i,yr]) ) ## trunate normal process error
      }
    }
    for(i in 2:q){ ##Done for the current year forecast. Excluded from pervious because n != q
      color[i,N] <- x[(i-1),N] + r * x[(i-1),N] * (1-x[(i-1),N])  ## latent process
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
    r ~ dexp(0.148) # Exp is the maximum entropy distribution for constraints of positive with given mean
    # 0.148 is from Richardson et al. 2006.
  }"
  }else{
    if(index=="GCC"){
      LogisticModel = "
      model{
      #### Data Models for complete years
      for(yr in 1:(N-1)){
      for(i in 1:n){
      p[i,yr] ~ dnorm(x[i,yr],p.PC)
      #mn[i,yr] ~ dnorm(x[i,yr],p.MN)
      #me[i,yr] ~ dnorm(x[i,yr],p.ME)
      }
      }
      ##Data Model for last year
      for(i in 1:q){
      p[i,N] ~ dnorm(x[i,N],p.PC)
      #mn[i,N] ~ dnorm(x[i,N],p.MN)
      #me[i,N] ~ dnorm(x[i,N],p.ME)
      }

      #### Process Model
      #### Color is the expected new phenology stage given the previous stage and logistic
      #### subtraction instead of addition in the discrete logistic eqn makes r negative (so logistic goes down).
      for(yr in 1:(N-1)){
      for(i in 2:n){
      color[i,yr] <- x[(i-1),yr] - r * x[(i-1),yr] * (1-x[(i-1),yr])  ## latent process
      xl[i,yr] ~ dnorm(color[i,yr],p.proc)  ## process error
      x[i,yr] <- max(0, min(1,xl[i,yr]) ) ## trunate normal process error
      }
      }
      for(i in 2:q){ ##Done for the current year forecast. Excluded from pervious because n != q
      color[i,N] <- x[(i-1),N] - r * x[(i-1),N] * (1-x[(i-1),N])  ## latent process
      xl[i,N] ~ dnorm(color[i,N],p.proc)  ## process error
      x[i,N] <- max(0, min(1,xl[i,N]) ) ## trunate normal process error
      }

      #### Priors
      for(yr in 1:N){ ##Initial Conditions
      x[1,yr] ~ dbeta(x1.a,x1.b)
      color[1,yr] ~ dbeta(x1.a,x1.b)
      }

      p.PC ~ dgamma(s1.PC,s2.PC)
      #p.MN ~ dgamma(s1.MN,s2.MN)
      #p.ME ~ dgamma(s1.ME,s2.ME)
      p.proc ~ dgamma(s1.proc,s2.proc)
      r ~ dexp(0.148) # Exp is the maximum entropy distribution for constraints of positive with given mean
      # 0.148 is from Richardson et al. 2006.
    }"
    }else{
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
      #### Color is the expected new phenology stage given the previous stage and logistic
      #### subtraction instead of addition in the discrete logistic eqn makes r negative (so logistic goes down).
      for(yr in 1:(N-1)){
        for(i in 2:n){
          color[i,yr] <- x[(i-1),yr] - r * x[(i-1),yr] * (1-x[(i-1),yr])  ## latent process
        xl[i,yr] ~ dnorm(color[i,yr],p.proc)  ## process error
        x[i,yr] <- max(0, min(1,xl[i,yr]) ) ## trunate normal process error
        }
      }
      for(i in 2:q){ ##Done for the current year forecast. Excluded from pervious because n != q
        color[i,N] <- x[(i-1),N] - r * x[(i-1),N] * (1-x[(i-1),N])  ## latent process
        xl[i,N] ~ dnorm(color[i,N],p.proc)  ## process error
        x[i,N] <- max(0, min(1,xl[i,N]) ) ## trunate normal process error
      }

      #### Priors
      for(yr in 1:N){ ##Initial Conditions
        x[1,yr] ~ dbeta(x1.a,x1.b)
        color[1,yr] ~ dbeta(x1.a,x1.b)
      }

      p.PC ~ dgamma(s1.PC,s2.PC)
      p.MN ~ dgamma(s1.MN,s2.MN)
      p.ME ~ dgamma(s1.ME,s2.ME)
      p.proc ~ dgamma(s1.proc,s2.proc)
      r ~ dexp(0.148) # Exp is the maximum entropy distribution for constraints of positive with given mean
      # 0.148 is from Richardson et al. 2006.
      }"
    }
  }


  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
