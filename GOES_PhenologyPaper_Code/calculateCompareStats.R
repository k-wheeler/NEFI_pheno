calTransRoots <- function(c,b,m,d){
  root1 <- (b*m+log(sqrt(3)+2))/b
  root2 <- (b*m+log(2-sqrt(3)))/b
  return(c(root1,root2))
}
calR2 <- function(obys,prys){
  obs.mean <- mean(obys)
  SStot <- sum((obys-obs.mean)**2)
  SSreg <- sum((prys-obs.mean)**2)
  SSres <- sum((obys-prys)**2)
  r2 <- (1-(SSres/SStot))
  if(r2<0){
    r2 <- 0
  }
  #adjr2 <- 1-(1-r2)*((length(obys)-1)/(length(obys)-2))
  return(r2)
}
calBiases <- function(xMatrix,yMatrix,name){
  bias <- matrix(nrow=0,ncol=3)
  allDif <- numeric()
  for(i in 1:nrow(xMatrix)){ ##each site
    dif <- mean((-xMatrix[i,]+yMatrix[i,]>0)) #changed signs
    allDif <- c(allDif,-xMatrix[i,]+yMatrix[i,]) #changed signs
    bias <- rbind(bias,quantile(-xMatrix[i,]+yMatrix[i,],c(0.025,0.5,0.975)))
  }
  
  ###bias<- rbind(bias,colMeans(bias)) ##<- this is wrong
  avg <- mean(allDif,na.rm = TRUE)
  avg.sd <- sqrt(var(allDif,na.rm=TRUE))
  avg.bot <- avg-1.96*(avg.sd/sqrt(length(allDif)))
  avg.top <- avg+1.96*(avg.sd/sqrt(length(allDif)))
  
  bias <- rbind(bias,c(avg.bot,avg,avg.top))
  rownames(bias) <- c(as.character(siteData$siteName[iseq]),"Average")
  fileName <- paste("PhenologyBias_",name,".csv",sep="")
  print(fileName)
  write.csv(file=fileName,bias,col.names = TRUE,row.names = TRUE,quote=F)
}

####Root Mean Square Error (RMSE)
calRMSE <- function(obs,preds){
  sums <- 0
  for(i in 1:length(obs)){
    sums <- sums+(preds[i]-obs[i])**2
  }
  return(sqrt(sums/length(obs)))
}

#' calculateCompareStates
#'
#' @param siteData 
#'
#' @return
#' @export
#'
#' @examples
calculateCompareStats <- function(siteData){
  n <- 10000
  ##1) Create data objects to be used in analyses and calculate transition dates
  startDay <- 1
  endDay <- 365
  
  iseq <- 1:nrow(siteData)
  yrs <- c(rep(2018,nrow(siteData)),rep(2019,nrow(siteData)-1))
  
  GOES.TranFs <- matrix(nrow=0,ncol=n)
  GOES.TranSs <- matrix(nrow=0,ncol=n)
  PC.TranFs <- matrix(nrow=0,ncol=n)
  PC.TranSs <- matrix(nrow=0,ncol=n)
  MODIS.N.TranFs <- matrix(nrow=0,ncol=n)
  MODIS.N.TranSs <- matrix(nrow=0,ncol=n)
  MODIS.E.TranFs <- matrix(nrow=0,ncol=n)
  MODIS.E.TranSs <- matrix(nrow=0,ncol=n)
  
  GOES.TranFs1 <- matrix(nrow=0,ncol=n)
  GOES.TranSs1 <- matrix(nrow=0,ncol=n)
  PC.TranFs1 <- matrix(nrow=0,ncol=n)
  PC.TranSs1 <- matrix(nrow=0,ncol=n)
  MODIS.N.TranFs1 <- matrix(nrow=0,ncol=n)
  MODIS.N.TranSs1 <- matrix(nrow=0,ncol=n)
  MODIS.E.TranFs1 <- matrix(nrow=0,ncol=n)
  MODIS.E.TranSs1 <- matrix(nrow=0,ncol=n)
  
  GOES.TranFs2 <- matrix(nrow=0,ncol=n)
  GOES.TranSs2 <- matrix(nrow=0,ncol=n)
  PC.TranFs2 <- matrix(nrow=0,ncol=n)
  PC.TranSs2 <- matrix(nrow=0,ncol=n)
  MODIS.N.TranFs2 <- matrix(nrow=0,ncol=n)
  MODIS.N.TranSs2 <- matrix(nrow=0,ncol=n)
  MODIS.E.TranFs2 <- matrix(nrow=0,ncol=n)
  MODIS.E.TranSs2 <- matrix(nrow=0,ncol=n)
  
  iVal=1
  
  for(iVal in 1:length(iseq)){
    i <- iseq[iVal]
    siteName <- as.character(siteData$siteName[i])
    print(siteName)
    yr <- yrs[iVal]
    URL <- as.character(siteData$URL[i])
    PFT <- as.character(siteData$PFT[i])
    lat <- as.character(siteData$Lat[i])
    long <- as.character(siteData$Long[i])
    TZ <- as.character(siteData$TZ[i])
    
    #GOES
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_GOES_varBurn.RData",sep="")
    load(inputFileName)
    var.mat <- as.matrix(var.burn)
    
    print(colnames(var.mat))
    GOES.TranF<-var.mat[,1]
    rndNums <- sample(1:length(GOES.TranF),n,replace=T)
    GOES.TranF <- GOES.TranF[rndNums]
    GOES.TranS<-var.mat[rndNums,2]
    GOES.bF <- var.mat[rndNums,3]
    GOES.bS <- var.mat[rndNums,4]
    GOES.c<-var.mat[rndNums,5]
    GOES.d<-var.mat[rndNums,6]
    GOES.TranFs <- rbind(GOES.TranFs,GOES.TranF)
    GOES.TranSs <- rbind(GOES.TranSs,GOES.TranS)
    
    #Calculate Roots (transition dates)
    rootsF1 <- numeric()
    rootsF2 <- numeric()
    rootsS1 <- numeric()
    rootsS2 <- numeric()
    for(g in 1:length(rndNums)){
      if(g %% 1000 == 0){
        print(g)
      }
      rootsF1<- c(rootsF1,calTransRoots(c=GOES.c[g],b=GOES.bF[g],m=GOES.TranF[g],d=GOES.d[g])[1])
      rootsF2<- c(rootsF2,calTransRoots(c=GOES.c[g],b=GOES.bF[g],m=GOES.TranF[g],d=GOES.d[g])[2])
      rootsS1<- c(rootsS1,calTransRoots(c=GOES.c[g],b=GOES.bS[g],m=GOES.TranS[g],d=GOES.d[g])[1])
      rootsS2<- c(rootsS2,calTransRoots(c=GOES.c[g],b=GOES.bS[g],m=GOES.TranS[g],d=GOES.d[g])[2])
    }
    GOES.TranFs1 <- rbind(GOES.TranFs1,rootsF2) #Naming switched for fall because of the formula. TranFs1 indicates the start of season
    GOES.TranFs2 <- rbind(GOES.TranFs2,rootsF1)
    GOES.TranSs1 <- rbind(GOES.TranSs1,rootsS1)
    GOES.TranSs2 <- rbind(GOES.TranSs2,rootsS2)
    
    #MODIS NDVI:
    if(yr==2018){
      startDate <- as.Date("2018-01-01")
      endDate <- as.Date("2018-12-31")
    }else if(yr==2019){
      startDate <- as.Date("2018-01-01")
      endDate <- as.Date("2018-12-31")
    }else{
      print("Error: year not recognized")
    }
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_NDVI_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.N.md.out)
    
    MODIS.N.TranF<-var.mat[,1]
    rndNums <- sample(1:length(MODIS.N.TranF),n,replace=T)
    MODIS.N.TranF <- MODIS.N.TranF[rndNums]
    MODIS.N.TranS<-var.mat[rndNums,2]
    MODIS.N.bF <- var.mat[rndNums,3]
    MODIS.N.bS <- var.mat[rndNums,4]
    MODIS.N.c<-var.mat[rndNums,5]
    MODIS.N.d<-var.mat[rndNums,6]
    MODIS.N.TranFs <- rbind(MODIS.N.TranFs,MODIS.N.TranF)
    MODIS.N.TranSs <- rbind(MODIS.N.TranSs,MODIS.N.TranS)
    
    rootsF1 <- numeric()
    rootsF2 <- numeric()
    rootsS1 <- numeric()
    rootsS2 <- numeric()
    for(g in 1:length(rndNums)){
      rootsF1<- c(rootsF1,calTransRoots(c=MODIS.N.c[g],b=MODIS.N.bF[g],m=MODIS.N.TranF[g],d=MODIS.N.d[g])[1])
      rootsF2<- c(rootsF2,calTransRoots(c=MODIS.N.c[g],b=MODIS.N.bF[g],m=MODIS.N.TranF[g],d=MODIS.N.d[g])[2])
      rootsS1<- c(rootsS1,calTransRoots(c=MODIS.N.c[g],b=MODIS.N.bS[g],m=MODIS.N.TranS[g],d=MODIS.N.d[g])[1])
      rootsS2<- c(rootsS2,calTransRoots(c=MODIS.N.c[g],b=MODIS.N.bS[g],m=MODIS.N.TranS[g],d=MODIS.N.d[g])[2])
    }
    MODIS.N.TranFs1 <- rbind(MODIS.N.TranFs1,rootsF2) #Naming switched for fall because of the formula. TranFs1 indicates the start of season
    MODIS.N.TranFs2 <- rbind(MODIS.N.TranFs2,rootsF1)
    MODIS.N.TranSs1 <- rbind(MODIS.N.TranSs1,rootsS1)
    MODIS.N.TranSs2 <- rbind(MODIS.N.TranSs2,rootsS2)
    
    print("MODIS NDVI Done")
    
    ##MODIS EVI
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_EVI_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(MODIS.E.md.out)
    
    MODIS.E.TranF<-var.mat[,1]
    rndNums <- sample(1:length(MODIS.E.TranF),n,replace=T)
    MODIS.E.TranF <- MODIS.E.TranF[rndNums]
    MODIS.E.TranS<-var.mat[rndNums,2]
    MODIS.E.bF <- var.mat[rndNums,3]
    MODIS.E.bS <- var.mat[rndNums,4]
    MODIS.E.c<-var.mat[rndNums,5]
    MODIS.E.d<-var.mat[rndNums,6]
    MODIS.E.TranFs <- rbind(MODIS.E.TranFs,MODIS.E.TranF)
    MODIS.E.TranSs <- rbind(MODIS.E.TranSs,MODIS.E.TranS)
    
    rootsF1 <- numeric()
    rootsF2 <- numeric()
    rootsS1 <- numeric()
    rootsS2 <- numeric()
    for(g in 1:length(rndNums)){
      rootsF1<- c(rootsF1,calTransRoots(c=MODIS.E.c[g],b=MODIS.E.bF[g],m=MODIS.E.TranF[g],d=MODIS.E.d[g])[1])
      rootsF2<- c(rootsF2,calTransRoots(c=MODIS.E.c[g],b=MODIS.E.bF[g],m=MODIS.E.TranF[g],d=MODIS.E.d[g])[2])
      rootsS1<- c(rootsS1,calTransRoots(c=MODIS.E.c[g],b=MODIS.E.bS[g],m=MODIS.E.TranS[g],d=MODIS.E.d[g])[1])
      rootsS2<- c(rootsS2,calTransRoots(c=MODIS.E.c[g],b=MODIS.E.bS[g],m=MODIS.E.TranS[g],d=MODIS.E.d[g])[2])
    }
    MODIS.E.TranFs1 <- rbind(MODIS.E.TranFs1,rootsF2) #Naming switched for fall because of the formula. TranFs1 indicates the start of season
    MODIS.E.TranFs2 <- rbind(MODIS.E.TranFs2,rootsF1)
    MODIS.E.TranSs1 <- rbind(MODIS.E.TranSs1,rootsS1)
    MODIS.E.TranSs2 <- rbind(MODIS.E.TranSs2,rootsS2)
    
    print("MODIS EVI Done")
    
    ##PC:
    inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_varBurn.RData",sep="")
    load(inputFileName)
    
    var.mat<-as.matrix(PC.md.out)
    
    PC.TranF<-var.mat[,1]
    rndNums <- sample(1:length(PC.TranF),n,replace=T)
    PC.TranF <- PC.TranF[rndNums]
    PC.TranS<-var.mat[rndNums,2]
    PC.bF <- var.mat[rndNums,3]
    PC.bS <- var.mat[rndNums,4]
    PC.c<-var.mat[rndNums,5]
    PC.d<-var.mat[rndNums,6]
    PC.TranFs <- rbind(PC.TranFs,PC.TranF)
    PC.TranSs <- rbind(PC.TranSs,PC.TranS)
    
    rootsF1 <- numeric()
    rootsF2 <- numeric()
    rootsS1 <- numeric()
    rootsS2 <- numeric()
    for(g in 1:length(rndNums)){
      rootsF1<- c(rootsF1,(calTransRoots(c=PC.c[g],b=PC.bF[g],m=PC.TranF[g],d=PC.d[g])[1]))
      rootsF2<- c(rootsF2,(calTransRoots(c=PC.c[g],b=PC.bF[g],m=PC.TranF[g],d=PC.d[g])[2]))
      rootsS1<- c(rootsS1,(calTransRoots(c=PC.c[g],b=PC.bS[g],m=PC.TranS[g],d=PC.d[g])[1]))
      rootsS2<- c(rootsS2,(calTransRoots(c=PC.c[g],b=PC.bS[g],m=PC.TranS[g],d=PC.d[g])[2]))
    }
    PC.TranFs1 <- rbind(PC.TranFs1,rootsF2) #Naming switched for fall because of the formula. TranFs1 indicates the start of season
    PC.TranFs2 <- rbind(PC.TranFs2,rootsF1)
    PC.TranSs1 <- rbind(PC.TranSs1,rootsS1)
    PC.TranSs2 <- rbind(PC.TranSs2,rootsS2)
    
    print("PC Done")
    
  }
  ##2) R2 and RMSE calculations
  IDs <- c("GOES_PC_spring50","MODISN_PC_spring50","MODISE_PC_spring50","MODISN_GOES_spring50","MODISE_GOES_spring50","MODISE_MODISN_spring50",
           "GOES_PC_fall50","MODISN_PC_fall50","MODISE_PC_fall50","MODISN_GOES_fall50","MODISE_GOES_fall50","MODISE_MODISN_fall50",
           "GOES_PC_springTran1","MODISN_PC_springTran1","MODISE_PC_springTran1","MODISN_GOES_springTran1","MODISE_GOES_springTran1","MODISE_MODISN_springTran1",
           "GOES_PC_fallTran1","MODISN_PC_fallTran1","MODISE_PC_fallTran1","MODISN_GOES_fallTran1","MODISE_GOES_fallTran1","MODISE_MODISN_fallTran1",
           "GOES_PC_springTran2","MODISN_PC_springTran2","MODISE_PC_springTran2","MODISN_GOES_springTran2","MODISE_GOES_springTran2","MODISE_MODISN_springTran2",
           "GOES_PC_fallTran2","MODISN_PC_fallTran2","MODISE_PC_fallTran2","MODISN_GOES_fallTran2","MODISE_GOES_fallTran2","MODISE_MODISN_fallTran2")
  
  R2s <- numeric()
  RMSEs <- numeric()
  xMatrix <- PC.TranSs
  yMatrix <- GOES.TranSs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranSs
  xMatrix <- PC.TranSs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.E.TranSs
  xMatrix <- PC.TranSs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranSs
  xMatrix <- MODIS.N.TranSs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranSs
  xMatrix <- MODIS.E.TranSs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranSs
  xMatrix <- MODIS.E.TranSs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  ##
  xMatrix <- PC.TranFs
  yMatrix <- GOES.TranFs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranFs
  xMatrix <- PC.TranFs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.E.TranFs
  xMatrix <- PC.TranFs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranFs
  xMatrix <- MODIS.N.TranFs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranFs
  xMatrix <- MODIS.E.TranFs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranFs
  xMatrix <- MODIS.E.TranFs
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  ######Tran 1
  xMatrix <- PC.TranSs1
  yMatrix <- GOES.TranSs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranSs1
  xMatrix <- PC.TranSs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.E.TranSs1
  xMatrix <- PC.TranSs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranSs1
  xMatrix <- MODIS.N.TranSs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranSs1
  xMatrix <- MODIS.E.TranSs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranSs1
  xMatrix <- MODIS.E.TranSs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  ##
  xMatrix <- PC.TranFs1
  yMatrix <- GOES.TranFs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranFs1
  xMatrix <- PC.TranFs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.E.TranFs1
  xMatrix <- PC.TranFs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranFs1
  xMatrix <- MODIS.N.TranFs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranFs1
  xMatrix <- MODIS.E.TranFs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranFs1
  xMatrix <- MODIS.E.TranFs1
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  #################Tran 2
  xMatrix <- PC.TranSs2
  yMatrix <- GOES.TranSs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranSs2
  xMatrix <- PC.TranSs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.E.TranSs2
  xMatrix <- PC.TranSs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranSs2
  xMatrix <- MODIS.N.TranSs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranSs2
  xMatrix <- MODIS.E.TranSs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranSs2
  xMatrix <- MODIS.E.TranSs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  ##
  xMatrix <- PC.TranFs2
  yMatrix <- GOES.TranFs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranFs2
  xMatrix <- PC.TranFs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.E.TranFs2
  xMatrix <- PC.TranFs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranFs2
  xMatrix <- MODIS.N.TranFs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- GOES.TranFs2
  xMatrix <- MODIS.E.TranFs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  yMatrix <- MODIS.N.TranFs2
  xMatrix <- MODIS.E.TranFs2
  R2s <- c(R2s,calR2(prys=rowMeans(xMatrix),obys=rowMeans(yMatrix)))
  RMSEs <- c(RMSEs,calRMSE(preds=rowMeans(yMatrix),obs=rowMeans(xMatrix)))
  
  output <- cbind(IDs,R2s,RMSEs)
  
  ##3) Bias
  ##Spring
  xMatrix <- PC.TranSs #x greater than y (Positive indicates that PC greater/later than GOES); Negative PC earlier than GOES; Negative GOES later than PC
  yMatrix <- GOES.TranSs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="GOES_PC_spring50_2")
  
  xMatrix <- PC.TranSs
  yMatrix <- MODIS.N.TranSs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_PC_spring50")
  
  xMatrix <- PC.TranSs
  yMatrix <- MODIS.E.TranSs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_PC_spring50")
  
  xMatrix <- GOES.TranSs
  yMatrix <- MODIS.E.TranSs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_GOES_spring50")
  
  xMatrix <- GOES.TranSs
  yMatrix <- MODIS.N.TranSs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_GOES_spring50")
  
  xMatrix <- MODIS.N.TranSs
  yMatrix <- MODIS.E.TranSs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_MODISN_spring50")
  
  xMatrix <- PC.TranSs1
  yMatrix <- GOES.TranSs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="GOES_PC_springTran1")
  
  xMatrix <- PC.TranSs1
  yMatrix <- MODIS.N.TranSs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_PC_springTran1")
  
  xMatrix <- PC.TranSs1
  yMatrix <- MODIS.E.TranSs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_PC_springTran1")
  
  yMatrix <- MODIS.E.TranSs1
  xMatrix <- GOES.TranSs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_GOES_springTran1")
  
  yMatrix <- MODIS.N.TranSs1
  xMatrix <- GOES.TranSs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_GOES_springTran1")
  
  xMatrix <- MODIS.N.TranSs1
  yMatrix <- MODIS.E.TranSs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_MODISN_springTran1")
  
  xMatrix <- PC.TranSs2
  yMatrix <- GOES.TranSs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="GOES_PC_springTran2")
  
  xMatrix <- PC.TranSs2
  yMatrix <- MODIS.N.TranSs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_PC_springTran2")
  
  xMatrix <- PC.TranSs2
  yMatrix <- MODIS.E.TranSs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_PC_springTran2")
  
  yMatrix <- MODIS.E.TranSs2
  xMatrix <- GOES.TranSs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_GOES_springTran2")
  
  yMatrix <- MODIS.N.TranSs2
  xMatrix <- GOES.TranSs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_GOES_springTran2")
  
  xMatrix <- MODIS.N.TranSs2
  yMatrix <- MODIS.E.TranSs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_MODISN_springTran2")
  
  ##Autumn
  xMatrix <- PC.TranFs
  yMatrix <- GOES.TranFs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="GOES_PC_fall50")
  
  xMatrix <- PC.TranFs
  yMatrix <- MODIS.N.TranFs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_PC_fall50")
  
  xMatrix <- PC.TranFs
  yMatrix <- MODIS.E.TranFs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_PC_fall50")
  
  yMatrix <- MODIS.E.TranFs
  xMatrix <- GOES.TranFs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_GOES_fall50")
  
  yMatrix <- MODIS.N.TranFs
  xMatrix <- GOES.TranFs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_GOES_fall50")
  
  xMatrix <- MODIS.N.TranFs
  yMatrix <- MODIS.E.TranFs
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_MODISN_fall50")
  
  
  xMatrix <- PC.TranFs1
  yMatrix <- GOES.TranFs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="GOES_PC_fallTran1")
  
  xMatrix <- PC.TranFs1
  yMatrix <- MODIS.N.TranFs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_PC_fallTran1")
  
  xMatrix <- PC.TranFs1
  yMatrix <- MODIS.E.TranFs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_PC_fallTran1")
  
  yMatrix <- MODIS.E.TranFs1
  xMatrix <- GOES.TranFs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_GOES_fallTran1")
  
  yMatrix <- MODIS.N.TranFs1
  xMatrix <- GOES.TranFs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_GOES_fallTran1")
  
  xMatrix <- MODIS.N.TranFs1
  yMatrix <- MODIS.E.TranFs1
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_MODISN_fallTran1")
  
  xMatrix <- PC.TranFs2
  yMatrix <- GOES.TranFs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="GOES_PC_fallTran2")
  
  xMatrix <- PC.TranFs2
  yMatrix <- MODIS.N.TranFs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_PC_fallTran2")
  
  xMatrix <- PC.TranFs2
  yMatrix <- MODIS.E.TranFs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_PC_fallTran2")
  
  yMatrix <- MODIS.E.TranFs2
  xMatrix <- GOES.TranFs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_GOES_fallTran2")
  
  yMatrix <- MODIS.N.TranFs2
  xMatrix <- GOES.TranFs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISN_GOES_fallTran2")
  
  xMatrix <- MODIS.N.TranFs2 ##x is second; y is first
  yMatrix <- MODIS.E.TranFs2
  calBiases(xMatrix = xMatrix,yMatrix = yMatrix, name="MODISE_MODISN_fallTran2")
  
  ##4) Load and add bias data to output table
  biasFiles <- dir(pattern="PhenologyBias")
  biases <- character()
  for(i in 1:length(IDs)){
    fileName <- paste("PhenologyBias_",IDs[i],".csv",sep="")
    dat <- read.csv(fileName,header=TRUE)

    newBias <- paste(round(dat[30,3],digits=2)," (",round((dat[30,4]-dat[30,2])/2,digits=2),")",sep="")
    biases <- c(biases,newBias)
  }
  output <- cbind(IDs,as.numeric(round(R2s,digits=2)),as.numeric(round(RMSEs,digits=2)),biases)
  colnames(output) <- c("ID","R^2","RMSE","Average Bias (95% CI)")
  
  write.table(output,file="GOES_Phenology_FitStatisticsUpdated.csv",row.names = FALSE,col.names = colnames(output),sep=",")
  
  ##5) Compare Transition Date Uncertainties
  GOES.CIs.TranSs1 <- cbind(data.frame(apply(GOES.TranSs1,1,quantile,0.975)-apply(GOES.TranSs1,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranSs1) <- c("CI","model")
  GOES.CIs.TranSs <- cbind(data.frame(apply(GOES.TranSs,1,quantile,0.975)-apply(GOES.TranSs,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranSs) <- c("CI","model")
  GOES.CIs.TranSs2 <- cbind(data.frame(apply(GOES.TranSs2,1,quantile,0.975)-apply(GOES.TranSs2,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranSs2) <- c("CI","model")
  
  PC.CIs.TranSs1 <- cbind(data.frame(apply(PC.TranSs1,1,quantile,0.975)-apply(PC.TranSs1,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranSs1) <- c("CI","model")
  PC.CIs.TranSs <- cbind(data.frame(apply(PC.TranSs,1,quantile,0.975)-apply(PC.TranSs,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranSs) <- c("CI","model")
  PC.CIs.TranSs2 <- cbind(data.frame(apply(PC.TranSs2,1,quantile,0.975)-apply(PC.TranSs2,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranSs2) <- c("CI","model")
  
  MODIS.E.CIs.TranSs1 <- cbind(data.frame(apply(MODIS.E.TranSs1,1,quantile,0.975)-apply(MODIS.E.TranSs1,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranSs1) <- c("CI","model")
  MODIS.E.CIs.TranSs <- cbind(data.frame(apply(MODIS.E.TranSs,1,quantile,0.975)-apply(MODIS.E.TranSs,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranSs) <- c("CI","model")
  MODIS.E.CIs.TranSs2 <- cbind(data.frame(apply(MODIS.E.TranSs2,1,quantile,0.975)-apply(MODIS.E.TranSs2,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranSs2) <- c("CI","model")
  
  MODIS.N.CIs.TranSs1 <- cbind(data.frame(apply(MODIS.N.TranSs1,1,quantile,0.975)-apply(MODIS.N.TranSs1,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranSs1) <- c("CI","model")
  MODIS.N.CIs.TranSs <- cbind(data.frame(apply(MODIS.N.TranSs,1,quantile,0.975)-apply(MODIS.N.TranSs,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranSs) <- c("CI","model")
  MODIS.N.CIs.TranSs2 <- cbind(data.frame(apply(MODIS.N.TranSs2,1,quantile,0.975)-apply(MODIS.N.TranSs2,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranSs2) <- c("CI","model")
  
  GOES.CIs.TranFs1 <- cbind(data.frame(apply(GOES.TranFs1,1,quantile,0.975)-apply(GOES.TranFs1,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranFs1) <- c("CI","model")
  GOES.CIs.TranFs <- cbind(data.frame(apply(GOES.TranFs,1,quantile,0.975)-apply(GOES.TranFs,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranFs) <- c("CI","model")
  GOES.CIs.TranFs2 <- cbind(data.frame(apply(GOES.TranFs2,1,quantile,0.975)-apply(GOES.TranFs2,1,quantile,0.025)),rep("GOES",29))
  colnames(GOES.CIs.TranFs2) <- c("CI","model")
  
  PC.CIs.TranFs1 <- cbind(data.frame(apply(PC.TranFs1,1,quantile,0.975)-apply(PC.TranFs1,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranFs1) <- c("CI","model")
  PC.CIs.TranFs <- cbind(data.frame(apply(PC.TranFs,1,quantile,0.975)-apply(PC.TranFs,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranFs) <- c("CI","model")
  PC.CIs.TranFs2 <- cbind(data.frame(apply(PC.TranFs2,1,quantile,0.975)-apply(PC.TranFs2,1,quantile,0.025)),rep("PC",29))
  colnames(PC.CIs.TranFs2) <- c("CI","model")
  
  MODIS.N.CIs.TranFs1 <- cbind(data.frame(apply(MODIS.N.TranFs1,1,quantile,0.975)-apply(MODIS.N.TranFs1,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranFs1) <- c("CI","model")
  MODIS.N.CIs.TranFs <- cbind(data.frame(apply(MODIS.N.TranFs,1,quantile,0.975)-apply(MODIS.N.TranFs,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranFs) <- c("CI","model")
  MODIS.N.CIs.TranFs2 <- cbind(data.frame(apply(MODIS.N.TranFs2,1,quantile,0.975)-apply(MODIS.N.TranFs2,1,quantile,0.025)),rep("MODISN",29))
  colnames(MODIS.N.CIs.TranFs2) <- c("CI","model")
  
  MODIS.E.CIs.TranFs1 <- cbind(data.frame(apply(MODIS.E.TranFs1,1,quantile,0.975)-apply(MODIS.E.TranFs1,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranFs1) <- c("CI","model")
  MODIS.E.CIs.TranFs <- cbind(data.frame(apply(MODIS.E.TranFs,1,quantile,0.975)-apply(MODIS.E.TranFs,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranFs) <- c("CI","model")
  MODIS.E.CIs.TranFs2 <- cbind(data.frame(apply(MODIS.E.TranFs2,1,quantile,0.975)-apply(MODIS.E.TranFs2,1,quantile,0.025)),rep("MODISE",29))
  colnames(MODIS.E.CIs.TranFs2) <- c("CI","model")
  
  CIs <- data.frame(rbind(PC.CIs.TranSs,GOES.CIs.TranSs,MODIS.N.CIs.TranSs,MODIS.E.CIs.TranSs))
  ###Checking to see if CI widths are significantly different than each other 
  pVals <- matrix(nrow=6,ncol=6)
  colnames(pVals) <- c("GOES and PC", "MODIS_N and PC","MODIS_E and PC",
                       "MODIS_N and GOES","MODIS_E and GOES","MODIS_E and MODIS_N")
  rowNum <- 1
  colNum <- 1
  
  cis1 <- PC.CIs.TranSs1[,1]
  cis2 <- GOES.CIs.TranSs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranSs1[,1]
  cis2 <- MODIS.N.CIs.TranSs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranSs1[,1]
  cis2 <- MODIS.E.CIs.TranSs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranSs1[,1]
  cis2 <- MODIS.N.CIs.TranSs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranSs1[,1]
  cis2 <- MODIS.E.CIs.TranSs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- MODIS.N.CIs.TranSs1[,1]
  cis2 <- MODIS.E.CIs.TranSs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  rowNum <- 2
  colNum <- 1
  
  cis1 <- PC.CIs.TranSs[,1]
  cis2 <- GOES.CIs.TranSs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranSs[,1]
  cis2 <- MODIS.N.CIs.TranSs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranSs[,1]
  cis2 <- MODIS.E.CIs.TranSs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranSs[,1]
  cis2 <- MODIS.N.CIs.TranSs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranSs[,1]
  cis2 <- MODIS.E.CIs.TranSs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- MODIS.N.CIs.TranSs[,1]
  cis2 <- MODIS.E.CIs.TranSs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  rowNum <- 3
  colNum <- 1
  
  cis1 <- PC.CIs.TranSs2[,1]
  cis2 <- GOES.CIs.TranSs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranSs2[,1]
  cis2 <- MODIS.N.CIs.TranSs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranSs2[,1]
  cis2 <- MODIS.E.CIs.TranSs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranSs2[,1]
  cis2 <- MODIS.N.CIs.TranSs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranSs2[,1]
  cis2 <- MODIS.E.CIs.TranSs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- MODIS.N.CIs.TranSs2[,1]
  cis2 <- MODIS.E.CIs.TranSs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  rowNum <- 4
  colNum <- 1
  
  cis1 <- PC.CIs.TranFs1[,1]
  cis2 <- GOES.CIs.TranFs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranFs1[,1]
  cis2 <- MODIS.N.CIs.TranFs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranFs1[,1]
  cis2 <- MODIS.E.CIs.TranFs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranFs1[,1]
  cis2 <- MODIS.N.CIs.TranFs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranFs1[,1]
  cis2 <- MODIS.E.CIs.TranFs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- MODIS.N.CIs.TranFs1[,1]
  cis2 <- MODIS.E.CIs.TranFs1[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  rowNum <- 5
  colNum <- 1
  
  cis1 <- PC.CIs.TranFs[,1]
  cis2 <- GOES.CIs.TranFs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranFs[,1]
  cis2 <- MODIS.N.CIs.TranFs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranFs[,1]
  cis2 <- MODIS.E.CIs.TranFs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranFs[,1]
  cis2 <- MODIS.N.CIs.TranFs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranFs[,1]
  cis2 <- MODIS.E.CIs.TranFs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- MODIS.N.CIs.TranFs[,1]
  cis2 <- MODIS.E.CIs.TranFs[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  rowNum <- 6
  colNum <- 1
  
  cis1 <- PC.CIs.TranFs2[,1]
  cis2 <- GOES.CIs.TranFs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranFs2[,1]
  cis2 <- MODIS.N.CIs.TranFs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- PC.CIs.TranFs2[,1]
  cis2 <- MODIS.E.CIs.TranFs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranFs2[,1]
  cis2 <- MODIS.N.CIs.TranFs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- GOES.CIs.TranFs2[,1]
  cis2 <- MODIS.E.CIs.TranFs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  colNum <- colNum + 1
  
  cis1 <- MODIS.N.CIs.TranFs2[,1]
  cis2 <- MODIS.E.CIs.TranFs2[,1]
  tTestVals <- t.test(cis1,cis2,paired = TRUE)
  pVals[rowNum,colNum] <- round(tTestVals$p.value,digits=3)
  
  write.table(pVals,"GOES_Phenology_CIwidth_pValuesUpdated.csv",sep=",",col.names = TRUE,row.names = FALSE)
}
