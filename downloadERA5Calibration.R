library("reticulate")
library(doParallel)
n.cores <- 16

#register the cores.
registerDoParallel(cores=n.cores)
setwd("/projectnb/dietzelab/kiwheel/ERA5")
siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)

cdsapi <- reticulate::import("cdsapi")
cclient <- cdsapi$Client()

# all_products <- c("reanalysis", "ensemble_members",
#                   "ensemble mean", "ensemble_spread")

end_date=as.Date("2019-12-31")

variables <- tibble::tribble(
  ~cf_name, ~units, ~api_name, ~ncdf_name,
  "air_temperature", "Kelvin", "2m_temperature", "t2m",
  "air_pressure", "Pa", "surface_pressure", NA_character_,
  NA_character_, "Kelvin", "2m_dewpoint_temperature", NA_character_,
  "precipitation_flux", "kg/m2/s", "total_precipitation", NA_character_,
  "eastward_wind", "m/s", "10m_u_component_of_wind", NA_character_,
  "northward_wind", "m/s", "10m_v_component_of_wind", NA_character_,
  "surface_downwelling_shortwave_flux_in_air", "W/m2", "surface_solar_radiation_downwards", NA_character_,
  "surface_downwelling_longwave_flux_in_air", "W/m2", "surface_thermal_radiation_downwards", NA_character_
)

var <- variables[["api_name"]][[1]]
foreach(i=2:nrow(siteData)) %dopar% {
#for(i in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  outfolder <- paste("Data/",siteName,sep="")
  dir.create(outfolder)
  lat <- as.numeric(siteData$Lat[i])
  long <- as.numeric(siteData$Long[i])
  start_date <- as.Date(siteData$startDate[i])
  
  area <- rep(round(c(lat, long) * 4) / 4, 2)
  ##Mean
  fname <- file.path(outfolder, paste(siteName,"_",start_date,"_",end_date,"_era5AirTemperatureMean.nc", sep =""))

  do_next <- tryCatch({
    cclient$retrieve(
      "reanalysis-era5-single-levels",
      list(
        variable = var,
        product_type = 'ensemble_mean',
        date = paste(start_date, end_date, sep = "/"),
        time = "00/to/23/by/1",
        area = area,
        grid = c(0.25, 0.25),
        format = "netcdf"
      ),
      fname
    )
    FALSE
  }, error = function(e) {
    print("Failed to download variable Mean")
    TRUE
  })
  
  fname <- file.path(outfolder, paste(siteName,"_",start_date,"_",end_date,"_era5AirTemperatureSpread.nc", sep =""))

  ##Spread
  do_next <- tryCatch({
    cclient$retrieve(
      "reanalysis-era5-single-levels",
      list(
        variable = var,
        product_type = 'ensemble_spread',
        date = paste(start_date, end_date, sep = "/"),
        time = "00/to/23/by/1",
        area = area,
        grid = c(0.25, 0.25),
        format = "netcdf"
      ),
      fname
    )
    FALSE
  }, error = function(e) {
    print("Failed to download variable Mean")
    TRUE
  })
}



# library(ncdf4)
# 
# tFile <- nc_open("harvard_2019-07-01_2019-07-03_era5AirTemperature.nc")
# Tairs <- ncvar_get(tFile)-273
# timeHours <- tFile$dim$time$vals #Hours since 1900-01-01 00:00:00.0 
# 
# ##Convert times to actual times 
# times <- as.POSIXct(timeHours*3600, origin = "1900-01-01",tz = "GMT")
# attributes(times)$tzone <- "America/New_York"
# plot(times,Tairs)