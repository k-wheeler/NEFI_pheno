library(PhenologyBayesModeling)
library(splines)
library(rjags)
library(PhenoForecast)
library(doParallel)
siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
#s <- 1
#detect cores.
n.cores <- 15

#register the cores.
registerDoParallel(cores=n.cores)

splineModel <- "
  model{
    ##Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(x[i],prec) ##Process error model
      x[i] <- muAvg + inprod(b[i,],beta[]) ##Process model

      yobs[i] ~ dnorm(y[i],obs.prec[i]) ##Data Model
    }
    ##Random effects
    for(i in 1:j){
      beta[i] ~ dnorm(0,taub)
    }

    ##Priors
    muAvg ~ dnorm(0,0.01)
    prec ~ dgamma(s1,s2)
    taub ~ dgamma(0.1,0.1)###0.1,0.1) ##??????? Not quite sure what this is yet
  }
"
nchain=3
startDate <- as.Date("2019-01-01")
endDate <- as.Date("2019-12-31")
output <- 
  foreach(s = 1:nrow(siteData)) %dopar% {
#for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[s])
  print(siteName)
  URL <- as.character(siteData$URL[s])
  
  fileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_EIV_varBurn_SF_spline.RData",sep="")
  if(!file.exists(fileName)){
    dat <- PC_data(siteName,URL,startDate,endDate,seasonOrder="SA",metric="GCC")
    datFinal <- list(yobs=dat$y,n=dat$n,obs.prec=dat$obs.prec)
    datFinal$j <- 30 #30 just looks so much better and gets the full increase of GCC at the beginning of summer
    datFinal$b <- bs(dat$x,datFinal$j)
    datFinal$s1 <- 0.001
    datFinal$s2 <- 0.00001
    
    j.model   <- jags.model(file = textConnection(splineModel),
                            data = datFinal,
                            n.chains = nchain)
    outBurn <- runForecastIter(j.model=j.model,variableNames = c("muAvg","prec","taub","x"),
                               baseNum = 100000,iterSize = 50000,effSize = 5000)
    ##Thin
    out.mat <- as.matrix(outBurn)
    thinAmount <- round(nrow(out.mat)/5000,digits=0)
    var.burn <- window(outBurn,thin=thinAmount)
    save(file=fileName,var.burn)
  }
}

# sum <- summary(outBurn20$predict)
# q20 <- sum$quantiles
# plot(dat$x,dat$y,pch=20)
# 
# lines(dat$x,q20[,1],col=2,lty=2)
# lines(dat$x,q20[,3],col=2,lty=1)
# lines(dat$x,q20[,5],col=2,lty=2)
# 
# lines(dat$x,q25[,1],col="blue",lty=2)
# lines(dat$x,q25[,3],col="blue",lty=1)
# lines(dat$x,q25[,5],col="blue",lty=2)
# 
# lines(dat$x,q30[,1],col="green",lty=2)
# lines(dat$x,q30[,3],col="green",lty=1)
# lines(dat$x,q30[,5],col="green",lty=2)
# 


