## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(GOESDiurnalNDVI)

## ------------------------------------------------------------------------
library(devtools)
#install_github("k-wheeler/NEFI_pheno/GOESDiurnalNDVI")

library(GOESDiurnalNDVI)


## ------------------------------------------------------------------------
###Site Characteristics:
siteName <- "RussellSage"
lat <- 32.457
long <- -91.9743
year <- 2017
TZ <- 5 #Time zone
day <- 234 #Day of year
dataPath <- "GOES_data" #Folder where the data is located
savePath <- paste(getwd(),"/",sep="")
print(savePath)
siteData <- cbind(siteName,lat,long,TZ)
fileName <- paste(savePath,"GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="")
if(!file.exists(fileName)){ #Note the file of the calculated NDVI's is included on github
  calculateNDVI_GOES_MAIN(day=day,siteData=siteData,year=year,TZ=TZ,
                          dataPath=dataPath,TZ_name="America/New_York",
                          savePath=savePath)
}

## ------------------------------------------------------------------------
fileName <- paste(savePath,"GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="")
GOESdat <- read.csv(fileName,header=FALSE)
plot(as.numeric(GOESdat[3,]),as.numeric(GOESdat[2,]),pch=20,xlim=c(0,20),ylim=c(0,1),
     ylab="NDVI",xlab="Time (Hour)")


## ------------------------------------------------------------------------
data <- list(x=as.numeric(GOESdat[3,]),y=as.numeric(GOESdat[2,]))
modelFitFileName <- paste(savePath,siteName,"_",year,day,"_varBurn.RData",sep="")
if(!file.exists(modelFitFileName)){
  j.model <- createDiurnalModel(siteName=siteName,data=data)
}


## ------------------------------------------------------------------------
if(!file.exists(modelFitFileName)){
  var.burn <- runMCMC_Model(j.model=j.model,variableNames=c("a","c","k","prec"),
                         baseNum=20000,iterSize =10000) #The baseNum and iterSize can be increased/decreased to make the code run faster if you know it will converge easier

save(var.burn,file=modelFitFileName)
}

## ------------------------------------------------------------------------
load(modelFitFileName)
plotCI(siteName=siteName,year=year,day=day,savePath=savePath)

## ------------------------------------------------------------------------
out.mat <- data.frame(as.matrix(var.burn))
c <- out.mat$c
c.quantiles <- quantile(c,c(0.025,0.50,0.975))
print(paste("Midday NDVI Quantiles (0.025,0.5,0.975):",as.character(c.quantiles[1]),
            as.character(c.quantiles[2]),as.character(c.quantiles[3])))
c.mean <- mean(c)
print(paste("Midday NDVI Mean Estimate:",as.character(c.mean)))

