##' Creates a logistic phenology forecast model based on PhenoCam and MODIS data
##'
##' @param data The data in the form of a list with data$p, data$mn, data$me, data$n, data$x_ic, and data$tau_ic
##' @param nchain The desired number of chains in the MCMC
##' @param season Season: spring or fall
##' @export
##' @import rjags
##' @import coda
logisticPhenoModel <- function(data,nchain,season){
  ##Set priors
  data$s1 <- 0.5
  data$s2 <- 0.2

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
  }
  }


  ###Create the JAGS model using the basic RandomWalk Model

  j.model   <- jags.model (file = textConnection(LogisticModel),
                           data = data,
                           n.chains = nchain)
  return(j.model)

}
