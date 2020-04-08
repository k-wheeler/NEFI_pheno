install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/GOESDiurnalNDVI",
                 repo=NULL,lib="/projectnb/dietzelab/kiwheel/Rlibrary")
library(GOESDiurnalNDVI)

savePath <- "PhenologyForecastData"
year <- 2019
dataPath <- paste("GOES_Data/GOES_Data",year,sep="") #Folder where the data is located
days <- seq(1,365)
mVersions <- c(rep(3,92),rep(6,273))
allSiteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)[1:19,]

TZ <- 5 #Time zone
siteData <- allSiteData[allSiteData$TZ==TZ,c(1,2,3,6)]
#i <- 1
for(i in 1:length(days)){
  mVersion <- mVersions[i]
  day <- days[i]
  print(day)
  calculateNDVI_GOES_MAIN(day=day,siteData=siteData,year=year,TZ=TZ,
                          dataPath=dataPath,TZ_name="America/New_York",
                          savePath=savePath,mVersion = mVersion)
}
