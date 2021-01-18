#' createAllDiurnalFits
#'
#' @return
#' @import PhenologyBayesModeling
#' @import GOESDiurnalNDVI
#' @import rjags
#' @import runjags
#' @import suncalc
#' @export
#'
#' @examples
createAllDiurnalFits <- function(calculatedNDVIGOESpath,DiurnalFitSavePath,siteData){
  # library("doParallel")
  # n.cores <- 10
  # registerDoParallel(cores=n.cores)
  days <- seq(1,365)
  for(year in c(2018,2019)){
    for(s in 1:nrow(siteData)){
      siteName <- as.character(siteData$siteName[s])
      lat <- as.numeric(siteData$Lat[s])
      long <- as.numeric(siteData$Long[s])
      TZ <- as.character(siteData$TZ[s])
      PFT <- as.character(siteData$PFT[s])
      print(siteName)    
      
      for(i in 1:length(days)){
        if(as.numeric(days[i])<10){
          days[i] <- paste("00",as.character(days[i]),sep="")
        }
        else if(as.numeric(days[i])<100){
          days[i] <- paste("0",days[i],sep="")
        }
      }
      
      #output <- foreach(day = days) %dopar% {
      for(day in days){
        fileName <- paste(calculatedNDVIGOESpath,"GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="") #Calculated NDVI file
        print(fileName)
        if(file.exists(fileName)){
          GOESdat <- read.csv(fileName,header=FALSE)
          # plot(as.numeric(GOESdat[3,]),as.numeric(GOESdat[2,]),pch=20,xlim=c(0,20),ylim=c(0,1),
          #      ylab="NDVI",xlab="Time (Hour)")
          if(ncol(GOESdat)>1){ #Need to check that it is not an empty file
            data <- list(x=as.numeric(GOESdat[3,]),y=as.numeric(GOESdat[2,]))
            modelFitFileName <- paste(DiurnalFitSavePath,siteName,"_",year,day,"_varBurn.RData",sep="")
            if(!file.exists(modelFitFileName)){
              j.model <- createDiurnalModel(siteName=siteName,data=data)
              var.burn <- runMCMC_Model(j.model=j.model,variableNames=c("a","c","k","prec"),
                                        baseNum=40000,iterSize =30000,maxGBR = 1.80) #The baseNum and iterSize can be increased/decreased to make the code run faster if you know it will converge easier
              ##Thin:
              if(typeof(var.burn)!=typeof(FALSE)){
                out.mat <- as.matrix(var.burn)
                thinAmount <- round(nrow(out.mat)/5000,digits=0)
                var.burn <- window(var.burn,thin=thinAmount)
                save(var.burn,file=modelFitFileName)
              }
            }
          }
        }
      }
    }
  }
}

