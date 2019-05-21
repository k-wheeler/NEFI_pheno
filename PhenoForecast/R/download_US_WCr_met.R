##' Downloads the meteological data for Willow Creek
##'
##' @param start_date The start date of desired data download
##' @param end_date The end date of desired data download
##' @import purrr
##' @export
download_US_WCr_met <- function(start_date, end_date) {
  print("Downloading Meterological Data")
  ##Following code provided by Katie Zarada
  #base_url <- "http://flux.aos.wisc.edu/data/cheas/wcreek/flux/prelim/wcreek" ##Old URL
  base_url <- "http://co2.aos.wisc.edu/data/cheas/wcreek/flux/prelim/"

  start_year <- lubridate::year(start_date)
  end_year <- lubridate::year(end_date)

  # Reading in the data
  raw.data <- start_year:end_year %>%
    purrr::map_df(function(syear) {
      influx <-
        tryCatch(
          read.table(
            paste0(base_url, syear, "_met.txt"),
            sep = "",
            header = FALSE
          ) %>%
            apply(2, trimws) %>%
            apply(2, as.character) %>%
            data.frame(stringsAsFactors = F),
          error = function(e) {
            NULL
          },
          warning = function(e) {
            NULL
          }
        )
    }) %>%
    mutate_all(funs(as.numeric))

  #Constructing the date based on the columns we have
  raw.data$date <-as.POSIXct(paste0(raw.data$V1,"/",raw.data$V2,"/",raw.data$V3," ", raw.data$V4 %>% as.integer(), ":",(raw.data$V4-as.integer(raw.data$V4))*60),
                             format="%Y/%m/%d %H:%M", tz="UTC")
  # Some cleaning and filtering
  raw.data <- raw.data %>%
    dplyr::select(V1,V2,V3,V4,V5, V6, V26, V35, V40, V59, date) %>%
    filter(date >= start_date & date <=end_date)

  #Colnames changed
  colnames(raw.data) <- c("Year", "Month", "Day", "Hour", "DoY", "FjDay", "Tair", "rH", "Tsoil", "Rg", "date")

  #####
  ##Subset and reformat to calculate average daily temperature
  dates <- as.Date(character())
  current.date <- as.Date(raw.data$date[1])
  Tairs <- numeric()
  sub.Tairs <- numeric()
  years <- numeric()
  months <- numeric()

  for(i in 1:nrow(raw.data)){
    if(as.Date(raw.data$date[i])==as.Date(current.date)){
      if(raw.data$Tair[i]!=-999){
        sub.Tairs <- c(sub.Tairs,raw.data$Tair[i])
      }
    }
    else{
      Tairs <- c(Tairs,mean(sub.Tairs))
      dates <- c(dates,current.date)
      months <- c(months,raw.data$Month[(i-1)])
      years <- c(years,raw.data$Year[(i-1)])
      current.date <- as.Date(raw.data$date[i])
      sub.Tairs <- numeric()
    }

  }
  #print(dates[1])
  #dat <- (cbind(dates,years,months,Tairs))
  dat <- data.frame(dates=dates,years=years,month=months,Tairs=Tairs)

  days <- seq(as.Date(start_date),as.Date(end_date),"day")
  print(length(days))
  newTairs <- rep(NA,length(days)) ##Averaged values
  for(i in 1:nrow(dat)){
    newTairs[which(days==dat[i,1])] <- dat[i,4]
  }
  newMonths <- lubridate::month(days)
  newYears <- lubridate::year(days)
  dat2 <- data.frame(dates=days,years=newYears,months=newMonths,Tairs=newTairs)

  #dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(1,181),]
  #dat2 <- dat2[dat2$months%in%seq(1,6,1),]
#
#   Tair <- matrix(nrow=181,ncol=0)
#   for(yr in lubridate::year(startDate):lubridate::year(endDate)){
#     subDat <- dat2[dat2$years==yr,]
#     print(length(subDat$Tairs))
#     Tair <- cbind(Tair,subDat$Tairs)
#   }
  return(dat2$Tairs)
}
