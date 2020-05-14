siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
savePath <- paste(getwd(),"/PhenologyForecastData/GOES_DiurnalFits/",sep="")
finalTotal <- 365+184
allFileCount <- 0
for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[s])
  print(siteName) 
  diurnal.files <- dir(path=savePath,pattern=siteName)
  print(length(diurnal.files))
  print(length(diurnal.files)/finalTotal)
  allFileCount <- allFileCount + length(diurnal.files)
}
print(allFileCount/(finalTotal*15))
  


# 
# siteName <- "howland"
# diurnal.files <- dir(pattern=paste("GOES_NDVI_Diurnal",siteName,"_",sep=""))
# 
# NDVI.files <- dir(pattern="GOES_NDVI_DiurnalrussellSage_")
# #length(NDVI.files)
# dat <- read.csv(NDVI.files[1],header=FALSE)
# for(i in 2:length(NDVI.files)){
#   print(NDVI.files[i])
#   dat <- cbind(dat,read.csv(NDVI.files[i],header=FALSE))
# }
# diurnal.files <- dir(pattern=paste("GOES_NDVI_Diurnal",siteName,"_",sep=""))
# days <- numeric()
# for(i in 1:length(diurnal.files)){
#   start <- as.numeric(strsplit(diurnal.files[i],"_")[[1]][4])
#   end <- as.numeric(strsplit(diurnal.files[i],"_")[[1]][5])
#   dys <- seq(start,end,1)
#   days <- c(days,dys)
# }
# #sort(days)
# all.days <- c(seq(1,321,1),seq(347,365,1))
# missingDays <- numeric()
# for(j in 1:length(all.days)){
#   if((!all.days[j] %in% days)){
#     missingDays <- c(missingDays,all.days[j])
#   }
# }
# missingDays
# length(missingDays)
