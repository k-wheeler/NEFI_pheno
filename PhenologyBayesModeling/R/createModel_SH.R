##' Create a Bayes Model for a shrubland site
##'
##' @param dataSource data source (GOES.NDVI, MODIS.NDVI, PC.GCC)
##' @param siteName Site Name
##' @param URL PhenoCam network URL
##' @param lat latitude of site in degrees
##' @param long longitude of site in degrees
##' @param startDay the day of year since 2017-01-01 to start the model
##' @param endDay the day of year since 2017-01-01 to end the model
##' @param niter the maximum number of iterations you want to give the model to converge within
##' @import rjags
##' @import runjags
##' @export
createBayesModel.SH <- function(dataSource,siteName="",URL="",startDay,endDay,lat,long,TZ=5) {
  nchain = 10
  inits <- list()
  if(dataSource=="PC.GCC"){
    data <- PC_data(siteName=siteName,URL=URL,startDay=startDay,endDay=endDay)
    inits.mu <- createInits(data=data,PFT="SH")
    if(siteName=="luckyHills"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,204,3),b=rnorm(1,-0.156,0.02),c=rnorm(1,0.07,0.02),d=rnorm(1,0.32,0.02),r=rnorm(1,-0.0156,0.002),k=rnorm(1,211,5))
        }
    }
    else if(siteName=="burns"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,150,10),b=rnorm(1,0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.02,0.002),k=rnorm(1,inits.mu$k,5))
      }
    }
    else{
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,220,3),b=rnorm(1,-0.13,0.02),c=rnorm(1,0.018,0.01),d=rnorm(1,0.3195,0.02),r=rnorm(1,-0.03,0.002),k=rnorm(1,240,5))
      }
    }
    print(inits)
    data$mean.d <- 0.3
    data$mean.c <- 0.1
  }
  else if(dataSource == "MODIS.NDVI"){
    data = MODIS_data(siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay,metric="NDVI")
    inits.mu <- createInits(data=data,PFT="SH")
    if(siteName=="luckyHills"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,204,3),b=rnorm(1,-0.156,0.02),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.0156,0.002),k=rnorm(1,211,5))
      }
    }
    else if(siteName=="burns"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,150,10),b=rnorm(1,-0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.02,0.002),k=rnorm(1,inits.mu$k,5))
      }
    }
    else{
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,220,3),b=rnorm(1,-0.13,0.02),c=rnorm(1,inits.mu$c,0.01),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.03,0.002),k=rnorm(1,240,5))
      }
    }
    data$mean.c <- 0.25
    data$mean.d <- 0.15
  }
  else if(dataSource == "MODIS.EVI"){
    data = MODIS_data(siteName=siteName,lat=lat,long=long,startDay = startDay,endDay = endDay,metric="EVI")
    inits.mu <- createInits(data=data,PFT="SH")
    if(siteName=="luckyHills"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,204,3),b=rnorm(1,-0.156,0.02),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.0156,0.002),k=rnorm(1,211,5))
      }
    }
    else if(siteName=="burns"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,150,10),b=rnorm(1,-0.11,0.05),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.02,0.002),k=rnorm(1,inits.mu$k,5))
      }
    }
    else{
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,220,3),b=rnorm(1,-0.13,0.02),c=rnorm(1,inits.mu$c,0.01),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.03,0.002),k=rnorm(1,240,5))
      }
    }
    data$mean.c <- 0.25
    data$mean.d <- 0.15
  }
  else if(dataSource=="GOES.NDVI"){
    data = GOES_data(siteName,startDay = startDay,endDay = endDay,lat=lat,long=long,TZ=TZ)
    inits.mu <- createInits(data=data,PFT="SH")
    if(siteName=="luckyHills"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,204,3),b=rnorm(1,-0.156,0.02),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.0156,0.002),k=rnorm(1,211,5))
      }
    }
    else if(siteName=="burns"){
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,150,1),b=rnorm(1,-0.14,0.01),c=rnorm(1,inits.mu$c,0.02),d=rnorm(1,inits.mu$d,0.001),r=rnorm(1,-0.02,0.002),k=rnorm(1,inits.mu$k,5))
      }
    }
    else{
      for(i in 1:nchain){
        inits[[i]] <- list(Tran=rnorm(1,220,3),b=rnorm(1,-0.13,0.02),c=rnorm(1,inits.mu$c,0.01),d=rnorm(1,inits.mu$d,0.02),r=rnorm(1,-0.03,0.002),k=rnorm(1,240,5))
      }
    }
    data$mean.c <- 0.25
    data$mean.d <- 0.15
  }

  data$mean.Tran <- 150
  data$p.Tran <- 1/(40**2)
  data$mean.b <- -0.10
  data$p.b <- 1/(0.05**2)
  data$p.c <- 1/(0.1**2)
  data$p.d <- 1/(0.1**2)
  data$mean.r <- -0.02
  data$p.r <- 1/(0.005**2)
  data$s1 <- 0.001
  data$s2 <- 0.00001
  data$mean.k <- 240
  data$p.k <- 1/(50**2)

  SH_model <- "
  model{
  ##priors
  Tran ~ dnorm(mean.Tran,p.Tran) ##Spring 50% transition date
  b ~ dnorm(mean.b,p.b) ##Fitting parameter for green up
  c ~ dnorm(mean.c,p.c) ##maximum-minimum
  d ~ dnorm(mean.d,p.d) ##minimum
  r ~ dnorm(mean.r,p.r) ##green down rate constant
  prec ~ dgamma(s1,s2)
  k ~ dnorm(mean.k,p.k) ##change point

  for(i in 1:n){
  mu1[i] <- c/(1+exp(b*(x[i]-Tran)))+d   	## process model for green up (logistic) (prebreak)
  mu2[i] <- c*exp(r*(x[i]-k))+d  ##process model for green down (exponential) (postbreak)
  mu[i] <- ifelse(x[i]>k,mu2[i],mu1[i])   #process model
  y[i]  ~ dnorm(mu[i],prec)		## data model
  }
  }
  "

  j.model   <- jags.model(file = textConnection(SH_model),
                          data = data,
                          inits=inits,
                          n.chains = nchain)
  return(j.model)
}
