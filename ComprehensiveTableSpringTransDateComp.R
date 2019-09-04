##Script to make the comprehensive plot/table comparing the predictability of the different models for different sites and different transitions
##Need to start with one site and one model and calculate how far off the transition date estimate occurs within the actual transition date CI

library(PhenoForecast)
library(PhenologyBayesModeling)
library(coda)
library(rjags)
library(doParallel)
library(ecoforecastR)
library("MODISTools")
library(doParallel)

##Set and register cores for parallel
n.cores <- 10
registerDoParallel(cores=n.cores)

#########Make phenofits for the last year of data
startDate <- as.Date("2018-07-01")
endDate <- as.Date("2019-06-30")
startDay <- 182
endDay <- 365+181
forecastLength <- 0
#iseq <- c(1,2,3,4,6,15,16,18,20,21,22,23,24)
#iseq <- c(1,2,3,4,6,15,16,20)
iseq <- c(1,2,3,4,6,15,16,18,20,24)
iseq <- c(1,2,3,4,6,15,16,20,24)
siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
forecastDataFolder <- "PhenologyForecastData/ForecastOutputs/AllForecasts/"
dataFolder <- "PhenologyForecastData/"
outTable <- matrix(nrow=0,ncol=6) ##Each col is a different model and rows are different sites
outTable2 <- matrix(nrow=0,ncol=3) ##Each col is a different model and rows are different sites
i <- 1
allDates <- c(seq(as.Date("2019-01-23"),as.Date("2019-01-25"),"day"),
              as.Date("2019-02-03"),as.Date("2019-02-05"),
              seq(as.Date("2019-02-07"),as.Date("2019-06-06"),"day"))
#output <- 
#  foreach(i=iseq) %dopar% {
for(i in iseq){  
  siteName <- as.character(siteData[i,1])
  print(siteName)
  siteStartDate <- as.character(siteData$startDate[i])
  outFileName <- paste("PhenologyForecastData/phenoFits/",siteName,"_forecast_spring2019_varBurn.RData",sep="")
  
  ##First need 50% transition date estimate from the model:
  load(outFileName)
  var.mat <- data.frame(as.matrix(varBurn))
  
  
  tran50 <- quantile(var.mat$TranS,c(0.025,0.5,0.975))-365
  plotDays <- seq(23,157)
  allPlotDays <- seq(-10,170)
  
  ##Random Walk
  outFileName <- paste(dataFolder, siteName, "_RW_spring50TranCompNEW.RData",sep="")
  if(!file.exists(outFileName)){
    allDates2 <- as.Date("2019-01-23")
    minsCI <- matrix(ncol=0,nrow=3) ##first column is 1-23, but the indexes start first of the year; 
      #Indices (values in minCI) are the day of the year; 22 needs to be added to the col number to get the day
    for(d in 1:length(allDates)){
      #for(d in 1:10){
      calEndDate <- allDates[d]
      #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
      #calEndDate <- as.Date("2019-05-20")
      print(calEndDate)
      forStartDate <- calEndDate + 1
      forEndDate <- (forStartDate+13)
      forDates <- seq(forStartDate,forEndDate,"day")
      yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
      
      forecastLength <- 14
      plotDates <- as.numeric(format(forDates,"%j"))
      
      RWfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_randomWalk_outBurn.RData",sep="")
      if(file.exists(RWfileName)){
        #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
        load(RWfileName)
        out.mat.par <- data.frame(as.matrix(outBurnRW$params))
        dayNumber <- dim(as.matrix(outBurnRW$predict))[2]-as.numeric(format(forEndDate,"%j"))+1
        out.mat <- as.matrix(outBurnRW$predict)
        out.mat <- out.mat[,dayNumber:ncol(out.mat)]
        mins <- numeric()
        for(s in 1:nrow(out.mat)){
          newVals <- which(out.mat[s,]>0.50)
          newMin <- min(newVals)
          mins <- c(mins,newMin)
        }
        withoutInfmins <- mins[!is.infinite(mins)]
        if(length(withoutInfmins)>(length(mins)/2)){
        
          minsCI <- cbind(minsCI,quantile(withoutInfmins,c(0.025,0.50,0.975)))
        }else{
          minsCI <- cbind(minsCI,c(NA,NA,NA))
        }
        allDates2 <- c(allDates2,calEndDate)
      }
    }
    allDates3 <- allDates2[2:length(allDates2)]
    RWdat <- list(dts=allDates3,minsCI=minsCI)
    save(RWdat,file = outFileName)
  }
  load(outFileName)
  allDates3 <- RWdat$dts
  minsCI <- RWdat$minsCI
  includeCol <- matrix(nrow=0,ncol=ncol(minsCI)) ##Columns to include (tran date close to overall tran CI)
  for(m in 1:ncol(minsCI)){
    includeCol[m] <- FALSE
    for(l in 1:3){
      if(!is.na(minsCI[l,m])){
      if(minsCI[l,m] < (tran50[3]+1) && minsCI[l,m] > (tran50[1]-1)){
        includeCol[m] <- TRUE
      }
      }
    }
  }
  potInd2 <- which(((minsCI[1,]!=Inf) * (minsCI[2,] != Inf) * (minsCI[3,]!=Inf))==1)
  potInd <- seq(1,ncol(minsCI))[includeCol] ##Indices for the columns that have a good tran date
  ##Need to see if the CI's overlap

  noFrstPot <- TRUE
  j <- 1
  frstDy <- min(potInd)
  # while(noFrstPot){ ##Only counts it as the first day if the next three days are also included
  #   if(potInd[j+1]>=potInd[j] &&
  #      potInd[j+2]>=potInd[j] &&
  #      potInd[j+3]>=potInd[j]){
  #     frstDy <- potInd[j] ##frstDy is the index in the potInd/minCIs
  #     noFrstPot <- FALSE
  #   }else{
  #     j <-  j + 1
  #   }
  #   if(j>(length(potInd)-3)){
  #     noFrstPot <- FALSE
  #     frstDy <- NA
  #   }
  # }

  noFrstPot <- TRUE
  j <- 1
  if(length(potInd2)>=3){
  while(noFrstPot){ ##Only counts it as the first day if the next three days are also included
    if(potInd2[j+1]>=potInd2[j] &&
       potInd2[j+2]>=potInd2[j] &&
       potInd2[j+3]>=potInd2[j]){
      frstDy2 <- potInd2[j]
      noFrstPot <- FALSE
    }else{
      j <-  j + 1
    }
    if(j>(length(potInd2)-3)){
      noFrstPot <- FALSE
      frstDy2 <- NA
    }
  }
  }else{
    frstDy2 <- min(potInd2)
  }




  newRow2 <- as.numeric(format(allDates3[frstDy2],"%j"))  ##Add 22 because the index counting starts 22 days late
  newRow <- as.numeric(format(allDates3[frstDy],"%j"))
  
  ###Logistic
  
  outFileName <- paste(dataFolder, siteName, "_L_spring50TranCompNEW.RData",sep="")
  if(!file.exists(outFileName)){
    allDates2 <- as.Date("2019-01-23")
    minsCI <- matrix(ncol=0,nrow=3)
    for(d in 1:length(allDates)){
      #for(d in 1:10){
      calEndDate <- allDates[d]
      #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
      #calEndDate <- as.Date("2019-05-20")
      print(calEndDate)
      forStartDate <- calEndDate + 1
      forEndDate <- (forStartDate+13)
      forDates <- seq(forStartDate,forEndDate,"day")
      yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
      
      forecastLength <- 14
      plotDates <- as.numeric(format(forDates,"%j"))
      LfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_logistic_outBurn.RData",sep="")
      if(file.exists(LfileName)){
        #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
        load(LfileName)
        out.mat.par <- data.frame(as.matrix(outBurnL$params))
        dayNumber <- dim(as.matrix(outBurnL$predict))[2]-as.numeric(format(forEndDate,"%j"))+1
        out.mat <- as.matrix(outBurnL$predict)
        out.mat <- out.mat[,dayNumber:ncol(out.mat)]
        mins <- numeric()
        for(s in 1:nrow(out.mat)){
          newVals <- which(out.mat[s,]>0.50)
          newMin <- min(newVals)
          mins <- c(mins,newMin)
        }
        withoutInfmins <- mins[!is.infinite(mins)]
        if(length(withoutInfmins)>(length(mins)/2)){
          
          minsCI <- cbind(minsCI,quantile(withoutInfmins,c(0.025,0.50,0.975)))
        }else{
          minsCI <- cbind(minsCI,c(NA,NA,NA))
        }
        allDates2 <- c(allDates2,calEndDate)
      }
    }
    allDates3 <- allDates2[2:length(allDates2)]
    Ldat <- list(dts=allDates3,minsCI=minsCI)
    save(Ldat,file = outFileName)
  }
  load(outFileName)

  allDates3 <- Ldat$dts
  minsCI <- Ldat$minsCI
  includeCol <- matrix(nrow=0,ncol=ncol(minsCI))
  for(m in 1:ncol(minsCI)){
    includeCol[m] <- FALSE
    
    # if(minsCI[1,m]!=Inf && minsCI[3,m]!= Inf){
    # ci <- seq(minsCI[1,m],minsCI[3,m],1)
    # }else if (minsCI[1,m]!=Inf && minsCI[2,m]!= Inf){
    #   ci <- seq(minsCI[1,m],minsCI[2,m],1)
    # }else if(minsCI[1,m]!=Inf){
    #   ci <- minsCI[1,m]
    # }else{ci <- FALSE}
    # if(typeof(ci)!=typeof(FALSE)){
    for(l in 1:3){
      if(!is.na(minsCI[l,m])){
      if(minsCI[l,m] < (tran50[3]+1) && minsCI[l,m] > (tran50[1]-1)){
        includeCol[m] <- TRUE
      }
    }
    }
    }
  potInd <- seq(1,ncol(minsCI))[includeCol]
  potInd2 <- which(((minsCI[1,]!=Inf) * (minsCI[2,] != Inf) * (minsCI[3,]!=Inf))==1)
  ##Need to see if the CI's overlap

  noFrstPot <- TRUE
  j <- 1
  frstDy <- min(potInd)
  # while(noFrstPot){ ##Only counts it as the first day if the next three days are also included
  #   if(potInd[j+1]>=potInd[j] &&
  #      potInd[j+2]>=potInd[j] &&
  #      potInd[j+3]>=potInd[j]){
  #     frstDy <- potInd[j]
  #     noFrstPot <- FALSE
  #   }else{
  #     j <-  j + 1
  #   }
  #   if(j>(length(potInd)-3)){
  #     noFrstPot <- FALSE
  #     frstDy <- NA
  #   }
  # }
  noFrstPot <- TRUE
  j <- 1
  if(length(potInd2)>=3){
  while(noFrstPot){ ##Only counts it as the first day if the next three days are also included
    if(potInd2[j+1]>=potInd2[j] &&
       potInd2[j+2]>=potInd2[j] &&
       potInd2[j+3]>=potInd2[j]){
      frstDy2 <- potInd2[j]
      noFrstPot <- FALSE
    }else{
      j <-  j + 1
    }
    if(j>(length(potInd2)-3)){
      noFrstPot <- FALSE
      frstDy2 <- NA
    }
  }
  }else{
    frstDy2 <- min(potInd2)
  }

  newRow2 <- c(newRow2,as.numeric(format(allDates3[frstDy2],"%j")))
  newRow <- c(newRow, as.numeric(format(allDates3[frstDy],"%j")))
  
  
  ###Logistic Covariate
  outFileName <- paste(dataFolder, siteName, "_LC_spring50TranCompNEW.RData",sep="")
  if(!file.exists(outFileName)){
    allDates2 <- as.Date("2019-01-23")
    minsCI <- matrix(ncol=0,nrow=3)
    for(d in 1:length(allDates)){
      #for(d in 1:10){
      calEndDate <- allDates[d]
      #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
      #calEndDate <- as.Date("2019-05-20")
      print(calEndDate)
      forStartDate <- calEndDate + 1
      forEndDate <- (forStartDate+13)
      forDates <- seq(forStartDate,forEndDate,"day")
      yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
      
      forecastLength <- 14
      plotDates <- as.numeric(format(forDates,"%j"))
      LCfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
      if(file.exists(LCfileName)){
        #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
        load(LCfileName)
        out.mat.par <- data.frame(as.matrix(outBurnLC2$params))
        dayNumber <- dim(as.matrix(outBurnLC2$predict))[2]-as.numeric(format(forEndDate,"%j"))+1
        out.mat <- as.matrix(outBurnLC2$predict)
        out.mat <- out.mat[,dayNumber:ncol(out.mat)]
        
        mins <- numeric()
        for(s in 1:nrow(out.mat)){
          newVals <- which(out.mat[s,]>0.50)
          newMin <- min(newVals)
          mins <- c(mins,newMin)
        }
        withoutInfmins <- mins[!is.infinite(mins)]
        if(length(withoutInfmins)>(length(mins)/2)){
          
          minsCI <- cbind(minsCI,quantile(withoutInfmins,c(0.025,0.50,0.975)))
        }else{
          minsCI <- cbind(minsCI,c(NA,NA,NA))
        }
        allDates2 <- c(allDates2,calEndDate)
      }
    }
    allDates3 <- allDates2[2:length(allDates2)]
    LCdat <- list(dts=allDates3,minsCI=minsCI)
    save(LCdat,file = outFileName)
    
  }
  load(outFileName)
  allDates3 <- LCdat$dts
  minsCI <- LCdat$minsCI

  includeCol <- matrix(nrow=0,ncol=ncol(minsCI))
  for(m in 1:ncol(minsCI)){
    includeCol[m] <- FALSE
    # if(minsCI[1,m]!=Inf && minsCI[3,m]!= Inf){
    #   ci <- seq(minsCI[1,m],minsCI[3,m],1)
    # }else if (minsCI[1,m]!=Inf && minsCI[2,m]!= Inf){
    #   ci <- seq(minsCI[1,m],minsCI[2,m],1)
    # }else if(minsCI[1,m]!=Inf){
    #   ci <- minsCI[1,m]
    # }else{ci <- FALSE}
    # if(typeof(ci)!=typeof(FALSE)){
      for(l in 1:3){
        if(!is.na(minsCI[l,m])){
        if(minsCI[l,m] < (tran50[3]+1) && minsCI[l,m] > (tran50[1]-1)){
          includeCol[m] <- TRUE
        }
        }
      }
    }
  potInd <- seq(1,ncol(minsCI))[includeCol]
  potInd2 <- which(((minsCI[1,]!=Inf) * (minsCI[2,] != Inf) * (minsCI[3,]!=Inf))==1)
  ##Need to see if the CI's overlap

  noFrstPot <- TRUE
  j <- 1
  frstDy <- min(potInd)
  # while(noFrstPot){ ##Only counts it as the first day if the next three days are also included
  #   if(potInd[j+1]>=potInd[j] &&
  #      potInd[j+2]>=potInd[j] &&
  #      potInd[j+3]>=potInd[j]){
  #     frstDy <- potInd[j]
  #     noFrstPot <- FALSE
  #   }else{
  #     j <-  j + 1
  #   }
  #   if(j>(length(potInd)-3)){
  #     noFrstPot <- FALSE
  #     frstDy <- NA
  #   }
  # }

  noFrstPot <- TRUE
  j <- 1
  if(length(potInd2)>=3){
  while(noFrstPot){ ##Only counts it as the first day if the next three days are also included
    if(potInd2[j+1]>=potInd2[j] &&
       potInd2[j+2]>=potInd2[j] &&
       potInd2[j+3]>=potInd2[j]){
      frstDy2 <- potInd2
      noFrstPot <- FALSE
    }else{
      j <-  j + 1
    }
    if(j>(length(potInd2)-3)){
      noFrstPot <- FALSE
      frstDy2 <- NA
    }
  }
  }else{
    frstDy2 <- min(potInd2)
  }

  newRow2 <- c(newRow2,as.numeric(format(allDates3[frstDy2],"%j")))
  newRow <- c(newRow, as.numeric(format(allDates3[frstDy],"%j")))

  newRow <- c(newRow,as.numeric(tran50))
  outTable <- rbind(outTable,newRow)
  outTable2 <- rbind(outTable2,newRow2)
}
outTable <- cbind(outTable,outTable[,5]-outTable[,1],
                  outTable[,5]-outTable[,2],
                  outTable[,5]-outTable[,3])

colnames(outTable) <- c("RW_DOY","L_DOY","LC_DOY","Actual_0.025","Actual_0.50",
                        "Actual_0.975","RW_diff","L_diff","LC_diff")
outTable <- data.frame(cbind(as.character(siteData$siteName[iseq]),outTable))
##Plotting
plotTable <- rbind(as.matrix(outTable$RW_diff),as.matrix(outTable$L_diff),as.matrix(outTable$LC_diff))
plotTable <- cbind(c(rep("A",length(iseq)),
                     rep("B",length(iseq)),
                     rep("C",length(iseq))),
                     plotTable)
colnames(plotTable) <- c("Model","Diff")

png(file="ESA_Presentation/springForecastHorizons.png", width=7, height=7,units="in",res=1000)
#par(mfrow=c(1,1),mai=c(1,1,1,0.5))
boxplot(as.numeric(plotTable[,2])~plotTable[,1],main="Forecast Horizon",lwd=2,cex.axis=2)
dev.off()


#########
newOutTable <- outTable[order(as.numeric(newOutTable$LC_diff)),]



