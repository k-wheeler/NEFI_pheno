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
    inits[[i]] <- list(a=rnorm(1,0.004,0.0001),c=rnorm(1,0.2,0.01))
  }
  data$alpha.c <- 2
  data$beta.c <- 1.5

  data$s1 <- 14
  data$s2 <- 0.13

  data$mean.a <- 0.0009
  data$p.a <- 1/(0.0003**2)
  data$alpha.p.cloud <- 1
  data$beta.p.cloud <- 17 #Set to have >90% less than 0.13 (based on a cloud mask accuracy of 0.87%)

  data$n <- length(data$x)
  ##Create k array

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
  #p.cloud ~ dunif(0,1)
  p.cloud ~ dbeta(alpha.p.cloud,beta.p.cloud)

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
