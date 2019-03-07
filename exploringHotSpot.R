i
library("suncalc")
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
i <- 9
siteName <- siteData$siteName[i]
lat <- siteData$Lat[i]
long <- siteData$Long[i]
#lat <- 0 #Bay of Fundy approximate
#long <- -75.2

##Calculating the satellite viewing angles

calSatAlt <- function(lat,long){
  ##Define constants
  satLat <- 0
  satLong <- -75.2
  
  R <- 6378140 #(in m from class notes)
  Rheight <- 35800*1000 #https://noaasis.noaa.gov/NOAASIS/ml/genlsatl.html
  Rsat <- R + Rheight
  
  term1 <- Rsat*(cos(satLong-long)*cos(satLat)*cos(lat)+sin(satLat)*sin(lat))-R
  term2 <- Rsat**2+R**2-2*R*Rsat*(cos(satLong - long)*cos(satLat)*cos(lat)+sin(satLat)*sin(lat))
  product <- term1 * term2**(-1/2)
  return(1.5708-acos(product))
}
calSatAzm <- function(lat,long){
  ##Define constants
  satLat <- 0
  satLong <- -75.2
  
  R <- 6378140 #(in m from class notes)
  Rheight <- 35800*1000 #https://noaasis.noaa.gov/NOAASIS/ml/genlsatl.html
  Rsat <- R + Rheight
  
  num <- sin(satLong-long)*cos(satLat)
  den <- cos(satLong - long)*cos(satLat)*sin(lat)-sin(satLat)*cos(lat)
  quot <- num/den
  return(atan(quot))
}
satAlt <- calSatAlt(lat=lat,long=long)
satAzm <- calSatAzm(lat=lat,long=long)

##Probably should loop through days 
#dates <- seq(as.Date("2017-07-01"),as.Date("2018-06-30"),"day")
#days <- as.character(seq(1,365,1))
dates <- seq(as.Date("2017-07-01"),as.Date("2017-07-31"),"day")
days <- format(dates,"%j")
for(d in 1:length(days)){
  if(as.numeric(days[d])<10){
    days[d] <- paste("00",days[d],sep="")
  }
  else if(as.numeric(days[d])<100){
    days[d] <- paste("0",days[d],sep="")
  }
}
alts <- numeric()
azms <- numeric()
timeVals <- numeric()
class(timeVals) <- "Date"

for(d in 1:length(days)){
  if(days[d]<182){
    yr <- 2018
  }
  else{
    yr <- 2017
  }
  yrDay <- paste("OR_ABI-L1b-RadC-M3C02_G16_s",yr,days[d],sep="")
  print(yrDay)
  files <- dir(path="GOES_Data2017",pattern=yrDay)
  for(f in 1:length(files)){
    hr <- substr(strsplit(files[f],split="_")[[1]][4],9,10)
    mn <- substr(strsplit(files[f],split="_")[[1]][4],11,12)
    timeVal <- paste(dates[f]," ",hr,":",mn,":00",sep="")
    pos <- getSunlightPosition(date=timeVal,lat=lat,lon=long)
    alts <- c(alts,pos$altitude)
    azms <- c(azms,pos$azimuth)
    timeVals <- c(timeVals,timeVal)
  }

}

outFileName <- paste(siteName,"_sunAngles.csv",sep="")
write.table(cbind(timeVals,alts,azms),file=outFileName,row.names=FALSE,col.names=TRUE)

# 
# ##Identify files
# days <- seq(as.Date("2018-01-01"),as.Date("2018-12-31"),"day")
# hrs <- c("00","01","02","03","04","05","06","07","08","09",as.character(seq(10,23)))
# mins <- c("00","01","02","03","04","05","06","07","08","09",as.character(seq(10,59)))
# 
# 
# 
# 
# alts2 <- numeric()
# azms2 <- numeric()
# 
# for(d in 1:length(days)){
#   print(d)
#   for(h in 1:length(hrs)){
#     #for(m in 1:length(mins)){
#     #timeVal <- paste(days[i]," ",hrs[h],":",mins[m],":00",sep="")
#       timeVal <- paste(days[i]," ",hrs[h],":00:00",sep="")
#       pos <- getSunlightPosition(date=timeVal,lat=lat,lon=long)
#       if(round(pos$altitude,digits=1)==round(satAlt,digits=1) && round(pos$azimuth,digits=1)==round(satAzm,digits=1)){
#         print(timeVal)
#         print(pos$altitude)
#         print(pos$azimuth)
#       }
#       
#       alts2 <- c(alts2,pos$altitude)
#       azms2 <- c(azms2,pos$azimuth)
#     #}
#   }
# }
# 
# #sunpos <- getSunlightPosition(date="2018-01-01 18:00:00",lat=lat,lon=long)
# 
# 
# print(range(azms2))
# print(satAzm)
# print(range(alts2))
# print(satAlt)
# which.min(abs(azms2-satAzm))
# which.min(abs(alts2-satAlt))
# which.min(abs(azms2-satAzm)+abs(alts2-satAlt))
# azms2[2]
# satAzm
# alts2[2]
# satAlt
