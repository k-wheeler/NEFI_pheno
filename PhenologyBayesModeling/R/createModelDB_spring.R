##' Create a spring Bayes Model for a deciduous broadleaf site
##'
##'@param data Data in the form of list that has data$y, data$x, and data$n defined
##' @param dataSource data source (GOES.NDVI, MODIS.NDVI, PC.GCC, MODIS.EVI)
##' @import rjags
##' @import runjags
##' @export
createModelDB_springFall <- function(data,dataSource) {
  nchain = 5
  if(dataSource =="PC.GCC"){
    data$mean.c <- 0.1
    data$mean.d <- 0.35
    data$p.c <- 1/(0.1**2)
    data$p.d <- 1/(0.1**2)
  }
  else if(dataSource == "MODIS.NDVI"){
    data$mean.c <- 0.4
    data$p.c <- 1/(0.2**2)
    data$mean.d <- 0.6
    data$p.d <- 1/(0.2**2)
  }
  else if(dataSource == "MODIS.EVI"){
    data$mean.c <- 0.4
    data$p.c <- 1/(0.2**2)
    data$mean.d <- 0.6
    data$p.d <- 1/(0.2**2)
  }
  else if(dataSource=="GOES.NDVI"){
    data$mean.c <- 0.4
    data$p.c <- 1/(0.2**2)
    data$mean.d <- 0.6
    data$p.d <- 1/(0.2**2)
  }
  data$s1 <- 0.001
  data$s2 <- 0.00001
  data$p.Tran <- 1/(40**2)
  data$p.b <- 1/(0.05**2)
  data$mean.TranS <- 475
  data$mean.bS <- -0.10
  data$mean.k <- 365
  data$p.k <- 1/(30**2)

  DB_model <- "
  model{
  ##priors
  TranS ~ dnorm(mean.TranS,p.Tran) ##S for spring
  bS ~ dnorm(mean.bS,p.b)
  d ~ dnorm(mean.d,p.d)
  c ~ dnorm(mean.c,p.c)
  k ~ dnorm(mean.k,p.k)
  prec ~ dgamma(s1,s2)

  for(i in 1:n){
  mu[i] <- c/(1+exp(bS*(x[i]-TranS)))+d ##process model for Spring

  y[i]  ~ dnorm(mu[i],prec)		## data model
  }
  }
  "
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          n.chains = nchain)

  return(j.model)
}
