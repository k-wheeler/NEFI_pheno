##' Calculates cummulative Tair (within one year for one ensemble)
##'
##' @param lat The site latitude
##' @param long The site longitude
##' @param years The desired years to download
##' @import xts
##' @export
load_ERA5_Tair <- function(lat,long,years) {
  source('/projectnb/dietzelab/hamzed/ERA5/nc_extractor_Ens.R')
  ##Could probably save this in a csv file so we don't have to process it all of the time
  out <- ERA5_extract_ENS(lat=lat, long=long,years=years,var="t2m")
  TairsOutput <- matrix(ncol=10,nrow=181*length(years))
  for(e in 1:10){
    t <- out[[e]]
    tDaily <- as.data.frame(xts::apply.daily(t,mean))
    tDaily <- tDaily %>% mutate(Date=rownames(.))
    #tSpring <- tDaily[(lubridate::month(tDaily[,2])%in% seq(1,6)),]
    tSpring <- tDaily[as.numeric(format(as.Date(tDaily[,2]),"%j"))%in% seq(1,181),]
    TairsOutput[,e] <- tSpring[,1]-273
  }
  return(TairsOutput)
}
