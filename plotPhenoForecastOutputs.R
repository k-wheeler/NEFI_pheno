install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenoForecast",repo=NULL)
library("PhenoForecast")
library("PhenologyBayesModeling")
library("coda")
library("dplyr")
library("rjags")


##Read in data
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="PhenologyForecastData/"
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"

forecastLength <- 15
#GEFS_files <- dir(path=dataDirectory,pattern="NOAA_GEFS")

endDate <- (Sys.Date()-1)
#startDate <- as.Date("2013-01-01")
#endDate <- as.Date("2019-01-27")
i <- 1
#i <- 1
siteName <- as.character(siteData[i,1])
print(siteName)
saveDirectory <- paste(dataDirectory,"ForecastOutputs/",siteName,"/",endDate,"/",sep="")
URL <- as.character(siteData[i,4])
lat <- as.numeric(siteData[i,2])
long <- as.numeric(siteData[i,3])
startDate <- as.Date(siteData[i,7])
outFileName <- paste("PhenologyForecastOutput_",siteName,"_",startDate,"_",endDate,".pdf",sep="")
pdf(outFileName,height=6,width=10)

##Load rescale data
rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_PC.csv",sep="")
rescaleData <- read.csv(rescaleFile,header=TRUE)
cMeans.p <- rescaleData$cMeans.p
dMeans.p <- rescaleData$dMeans.p
rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_MN.csv",sep="")
rescaleData <- read.csv(rescaleFile,header=TRUE)
cMeans.mn <- rescaleData$cMeans.mn
dMeans.mn <- rescaleData$dMeans.mn
rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_ME.csv",sep="")
rescaleData <- read.csv(rescaleFile,header=TRUE)
cMeans.me <- rescaleData$cMeans.me
dMeans.me <- rescaleData$dMeans.me

##Random Walk:
RWFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_randomWalk_outBurn.RData",sep="")
load(RWFile)

##Plot Parameter Density Plots
out.mat.par <- data.frame(as.matrix(outBurnRW$params))
colnames(out.mat.par)
par(mfrow=c(2,2))
plot(density(out.mat.par$p.ME),main="Density of p.ME")
plot(density(out.mat.par$p.MN),main="Density of p.MN")
plot(density(out.mat.par$p.PC),main="Density of p.PC")
plot(density(out.mat.par$p.proc),main="Density of p.proc")

plot(density(sqrt(1/out.mat.par$p.PC)),main="Density of PC Obs. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.ME)),main="Density of ME Obs. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.MN)),main="Density of MN Ob. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.proc)),main="Density of Process Error (SD)")

par(mfrow=c(1,1))
##Plot Forecast including calibration
dayNumber <- dim(as.matrix(outBurnRW$predict))[2]
out.mat.RW <- as.matrix(outBurnRW$predict)
plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.RW,forecastType = "randomWalk",days=seq(1,dayNumber,1))
offset <- 365-as.numeric(format(startDate,"%j"))
for(i in seq(offset,dayNumber,365)){
  abline(v=i,col="red")
}
abline(v=(dayNumber-forecastLength+1),col="purple")


##Plot Current year
lengthLastYear <- (as.numeric(format(endDate,"%j"))+forecastLength)
print(lengthLastYear)
lastYearIndices <- seq(((dayNumber-lengthLastYear)+2),dayNumber,1)
out.mat.lastYear <- out.mat.RW[,lastYearIndices]

plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.lastYear,forecastType = "randomWalk",days=seq(1,lengthLastYear,1))
abline(v=(lengthLastYear-forecastLength+1),col="purple")

##Add on data:
##PhenoCam
#PCfileName <- paste(dataDirectory,siteName,"_",startDate,"_",endDate,"_PC_Data.RData",sep="")
#print(PCfileName)
#print(file.exists(PCfileName))
#if(!file.exists(PCfileName)){
  print(URL)
  phenoData <- download.phenocam(URL)
  print(dim(phenoData))
  #save(phenoData,file=PCfileName)
#}
#load(PCfileName)
p <- phenoData$gcc_mean[phenoData$year==2019]

print(p)
time.p.date <- as.Date(phenoData$date[phenoData$year==2019])
time.p <-  as.numeric(format(time.p.date,"%j"))

p <- rescale(c=cMeans.p[length(cMeans.p)],d=dMeans.p[length(dMeans.p)],yseq=p)

print(length(time.p))
print(length(p))
points(time.p,p,pch=20,col="red")

##MODIS NDVI and EVI
mn <- prepareMODIS(startDate=startDate,endDate=endDate,metric="NDVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
me <- prepareMODIS(startDate=startDate,endDate=endDate,metric="EVI",timeForecast=time.p.date,dataDirectory=dataDirectory,siteName=siteName)
points(time.p,mn,col="green",pch="+")
points(time.p,me,col="green",pch=1)
legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))

##Basic logistic:
LFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_logistic_outBurn.RData",sep="")
load(LFile)
out.mat.par <- data.frame(as.matrix(outBurnL$params))
colnames(out.mat.par)
par(mfrow=c(2,3))
plot(density(out.mat.par$p.ME),main="Density of p.ME")
plot(density(out.mat.par$p.MN),main="Density of p.MN")
plot(density(out.mat.par$p.PC),main="Density of p.PC")
plot(density(out.mat.par$p.proc),main="Density of p.proc")
plot(density(out.mat.par$r),main="Density of r")
par(mfrow=c(2,3))
plot(density(sqrt(1/out.mat.par$p.PC)),main="Density of PC Obs. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.ME)),main="Density of ME Obs. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.MN)),main="Density of MN Ob. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.proc)),main="Density of Process Error (SD)")
plot(density(out.mat.par$r),main="Density of r")

##Plot
par(mfrow=c(1,1))
dayNumber <- dim(as.matrix(outBurnL$predict))[2]
out.mat.L <- as.matrix(outBurnL$predict)
plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.L,forecastType = "Logistic",days=seq(1,dayNumber,1))#2484

for(i in seq(1,2484,182)){
  abline(v=i,col="red")
}
abline(v=(dayNumber-forecastLength+1),col="purple")
#points(time.p,p,pch=20,col="red")
#points(time.p,mn,col="green",pch=3)
#points(time.p,me,col="green",pch=1)
legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))


##Plot Current year
lastYearIndices <- seq(((dayNumber-lengthLastYear)+2),dayNumber,1)
out.mat.lastYear <- out.mat.L[,lastYearIndices]

plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.lastYear,forecastType = "Logistic",days=seq(1,lengthLastYear,1))
abline(v=(lengthLastYear-forecastLength+1),col="purple")

##Add on data points
points(time.p,p,pch=20,col="red")
points(time.p,mn,col="green",pch=3)
points(time.p,me,col="green",pch=1)
legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))

##Logistic with covariates model
LCFile <- paste(saveDirectory,siteName,"_",startDate,"_",endDate,"_LC_outBurn.RData",sep="")
load(LCFile)
out.mat.par <- data.frame(as.matrix(outBurnLC$params))
colnames(out.mat.par)
par(mfrow=c(2,3))
plot(density(out.mat.par$p.ME),main="Density of p.ME")
plot(density(out.mat.par$p.MN),main="Density of p.MN")
plot(density(out.mat.par$p.PC),main="Density of p.PC")
plot(density(out.mat.par$p.proc),main="Density of p.proc")
plot(density(out.mat.par$b1),main="Density of b1")
plot(density(out.mat.par$b0),main="Density of b0")

par(mfrow=c(2,3))
plot(density(sqrt(1/out.mat.par$p.PC)),main="Density of PC Obs. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.ME)),main="Density of ME Obs. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.MN)),main="Density of MN Ob. Error (SD)")
plot(density(sqrt(1/out.mat.par$p.proc)),main="Density of Process Error (SD)")

plot(density(out.mat.par$b1),main="Density of b1")
plot(density(out.mat.par$b0),main="Density of b0")

##Plot
par(mfrow=c(1,1))
dayNumber <- dim(as.matrix(outBurnLC$predict))[2]
out.mat.LC <- as.matrix(outBurnLC$predict)
plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.LC,forecastType = "Logistic Covariate",days=seq(1,dayNumber,1))#2484

for(i in seq(1,2484,182)){
  abline(v=i,col="red")
}
abline(v=(dayNumber-forecastLength+1),col="purple")

##Plot Current year
lastYearIndices <- seq(((dayNumber-lengthLastYear)+2),dayNumber,1)
out.mat.lastYear <- out.mat.LC[,lastYearIndices]

plotForecastOutput(siteName=siteName,URL=URL,forecastLength=forecastLength,out.mat=out.mat.lastYear,forecastType = "Logistic Covariate",days=seq(1,lengthLastYear,1))
abline(v=(lengthLastYear-forecastLength+1),col="purple")

##Add on data points
points(time.p,p,pch=20,col="red")
points(time.p,mn,col="green",pch=3)
points(time.p,me,col="green",pch=1)
legend("topleft",c("PC","MODIS NDVI","MODIS EVI"),col=c("red","green","green"),pch=c(20,3,1))

##Plot Sf vs R




dev.off()
