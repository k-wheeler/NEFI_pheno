##' Create a diurnal model for a deciduous broadleaf site
##'
##' @param siteName Site Name
##' @param data Data object
##' @import rjags
##' @import runjags
##' @export
createBayesModel.Diurnal <- function(siteName,data){
  nchain <-  5
  inits <- list()

  for(i in 1:5){
    inits[[i]] <- list(a=rnorm(1,0.004,0.0001),c=rnorm(1,0.2,0.01),k=rnorm(1,12,0.2))
  }
  #data$y <- data$y+0.2
  print("Updated")
  #data$mean.c <- 0.75
  #data$p.c <- 1/(0.05**2)
  data$alpha.c <- 2
  data$beta.c <- 1.5
  #data$min.c <- 0
  #data$max.c <- 1
  #data$s1 <- 0.001 #0.7
  #data$s2 <- 0.00001 #0.08
  data$s1 <- 14
  data$s2 <- 0.13
  #data$p.Tran <- 1/(1**2)
  #data$p.b <- 1/(1**2)
  #data$mean.TranL <- 0# 7.5
  #data$mean.bL <- -1.5
  #data$mean.TranR <- 25#17.5
  #data$mean.bR <- 1.8
  data$mean.a <- 0.0009
  data$p.a <- 0.0003
  data$mean.k <- 12
  data$p.k <- 1/(1**2)
  data$n <- length(data$x)
  print("finished defining data")
  DB_model_MM <- "
  model{
  ##priors
  #TranL ~ dnorm(mean.TranL,p.Tran) ##S for spring
  #bL ~ dnorm(mean.bL,p.b)
  #TranR ~ dnorm(mean.TranR,p.Tran)  ##F for fall/autumn
  #bR ~ dnorm(mean.bR,p.b)
  a ~ dnorm(mean.a,p.a) I(0,)
  c ~ dbeta(alpha.c,beta.c)
  #c ~ dunif(min.c,max.c)
  k ~ dnorm(mean.k,p.k)
  prec ~ dgamma(s1,s2)
  alp ~ dunif(1,100)
  bet ~ dunif(1,100)
  p.cloud ~ dunif(0,1)

  for(i in 1:n){
  muL[i] <- -a * exp(-1 * (x[i]-k)) + c + a
  muR[i] <- -a * exp((x[i]-k)) + c + a

  f[i] <- ifelse(x[i]>k,muR[i],muL[i])   #change point process model

  y[i] ~ dnorm(mu[i],prec)   ##data model
  is.cloudy[i] ~ dbern(p.cloud)
  trans[i] ~ dbeta(alp,bet)
  mu[i] <- is.cloudy[i] * trans[i]*f[i] + (1-is.cloudy[i]) * f[i]

  }
  }
  "

  j.model   <- jags.model(file = textConnection(DB_model_MM),
                          data = data,
                          inits = inits,
                          n.chains=nchain)
  return(j.model)

}
