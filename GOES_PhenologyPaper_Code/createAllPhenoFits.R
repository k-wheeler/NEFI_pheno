#' createAllPhenoFits
#'
#' @param startDate 
#' @param endDate 
#' @param siteData 
#'
#' @return
#' @import PhenologyBayesModeling
#' @import rjags
#' @import runjags
#' @import MODISTools
#' @import ncdf4
#' @import plyr
#' @export
#'
#' @examples
createAllPhenoFits <- function(startDate,endDate,siteData){

  DB.vars <- c("TranF","bF","TranS","bS","c","d","prec")
  # library(doParallel)
  # n.cores <- 16
  # registerDoParallel(cores=n.cores)
  # output <- 
  #   foreach(i=1:nrow(siteData)) %dopar% {
    for(i in 1:nrow(siteData)){
      siteName <- as.character(siteData$siteName[i])
      print(siteName)
      URL <- as.character(siteData$URL[i])
      PFT <- as.character(siteData$PFT[i])
      print(PFT)
      lat <- as.character(siteData$Lat[i])
      long <- as.character(siteData$Long[i])
      TZ <- as.character(siteData$TZ[i])
      ##PhenoCam:
      #fileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_EIV_varBurn_SF.RData",sep="")
      fileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_varBurn.RData",sep="")
      if(!file.exists(fileName)){
        j.model.PC <- createBayesModel.DB(dataSource="PC.GCC",siteName=siteName,URL=URL,seasonOrder = "SF",startDate = startDate,
                                          endDate=endDate)
        PC.md.out <- runMCMC_Model(j.model=j.model.PC,variableNames = DB.vars,baseNum = 5000,iterSize = 5000,maxGBR = 5)
        save(PC.md.out,file=fileName)
      }
      
      ##MODIS NDVI
      metric <- "NDVI"
      #fileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_NDVI_EIV_varBurn.RData",sep="")
      fileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_NDVI_varBurn.RData",sep="")
      print(fileName)
      if(!file.exists(fileName)){
        j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.NDVI",siteName=siteName,lat=lat,long=long,
                                             startDate = startDate,endDate = endDate,seasonOrder = "SF")
        MODIS.N.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars,baseNum=5000,iterSize = 5000)
        save(MODIS.N.md.out,file=fileName)
      }
      
      ##MODIS EVI
      #fileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_EVI_EIV_varBurn.RData",sep="")
      fileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_EVI_varBurn.RData",sep="")
      print(fileName)
      metric <- "EVI"
      if(!file.exists(fileName)){
        print("MODIS File Downloaded")
        j.model.MODIS <- createBayesModel.DB(dataSource="MODIS.EVI",siteName=siteName,lat=lat,long=long,
                                             startDate = startDate,endDate = endDate,seasonOrder = "SF")
        MODIS.E.md.out <- runMCMC_Model(j.model=j.model.MODIS,variableNames = DB.vars,baseNum=5000,iterSize = 5000)
        save(MODIS.E.md.out,file=fileName)
      }
      
      ##GOES NDVI
      varBurnFileName <- paste(siteName,"_",startDate,"_",endDate,"_GOES_varBurn.RData",sep="")
      if(!file.exists(varBurnFileName)){
        outDataFile <- paste(siteName,"_",yr,"_diurnalFitDataFiltered.RData",sep="") 
        load(outDataFile)
        j.model <- createBayesModel.DB_Overall(data=filteredData,seasonOrder="SF")
        var.burn <- runMCMC_Model(j.model=j.model,variableNames = DB.vars,baseNum=50000,iterSize=20000)
        save(var.burn,file=varBurnFileName)
      }
    }
  
}