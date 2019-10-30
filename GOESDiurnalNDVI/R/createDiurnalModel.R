#' Creates the diurnal jags model object
#'
#' @param siteName The site name used for file naming
#' @param data List that includes the data to be fitted (e.g. data$x and data$y) where x is the vector of hours of the day (minutes and seconds included in the decimal) y is the vector of NDVI values
#' @export
#' @import rjags
#' @import runjags
#' @import coda
createDiurnalModel <- function(siteName,data){
  print("entered model")
  nchain <-  5
  data$alpha.c <- 2
  data$beta.c <- 1.5
  data$s1 <- 0.001
  data$s2 <- 0.00001

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
                          n.chains=nchain)
  return(j.model)

}
