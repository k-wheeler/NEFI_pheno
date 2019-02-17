library("rjags")
library("runjags")

##' Create a Bayes Model for an evergreen needleleaf site
##'
##' @param dataSource data source (GOES.NDVI, MODIS.NDVI, PC.GCC)
##' @param siteName Site Name
##' @param URL PhenoCam network URL
##' @import rjags
##' @import runjags
##' @export
createBayesModel.EN <- function(dataSource,siteName="",URL="") {
  ##NOTE: Haven't gotten this to work well
  #siteName <- "niwot"
  dataSource <- "PC.GCC"
  nchain = 10
  if(dataSource=="PC.GCC"){
    fileName <- paste(siteName,"_PC_Data.RData",sep="")
    load(fileName)
    data <- PC.data
    #data = PC_data(URL=URL,startDay = 110,endDay = 474)
    data$mean.a <- 0.36
    data$p.a <- 0.05
    data$mean.d <- 0.33
    data$p.d <- 0.05
  }
  else if(dataSource == "MODIS.NDVI"){
    data = MODIS_data(siteName=siteName)
    data$mean.a <- 0.6
    data$p.a <- 1/(0.1**2)
    data$mean.d <- 0.3
    data$p.d <- 1/(0.1**2)
  }
  else if(dataSource=="GOES.NDVI"){
    data = GOES_data(siteName,startDay = 110,endDay = 474)
    data$mean.a <- 0.6
    data$p.a <- 1/(0.1**2)
    data$mean.d <- 0.3
    data$p.d <- 1/(0.1**2)
  }
  inits <- list()
  inits.mu <- createInits_EN(data)
  for(i in 1:nchain){
    inits[[i]] <- list(a=rnorm(1,inits.mu$a,0.02),b=rnorm(1,inits.mu$b,4),d=rnorm(1,inits.mu$d,0.02))
  }

  data$mean.b <- 75
  data$p.b <- 3
  data$s1.c <- 0.001
  data$s2.c <- 0.00001
  data$s1.prec <- 0.001
  data$s2.prec <- 0.00001

  EN_model <- "
  model{
  ##priors
  a ~ dnorm(mean.a,p.a)
  b ~ dnorm(mean.b,p.b)
  c ~ dgamma(s1.c,s2.c)
  d ~ dnorm(mean.d,p.d)
  prec ~ dgamma(s1.prec,s2.prec)

  for(i in 1:n){
  mu[i] <- a*exp((-1*(x[i]-b)**2)/(2*c**2))+d   	## process model
  y[i]  ~ dnorm(mu[i],prec)
  }
  }
  "

  j.model   <- jags.model(file = textConnection(EN_model),
                          data = data,
                          n.chains = nchain)

  return(j.model)
}
