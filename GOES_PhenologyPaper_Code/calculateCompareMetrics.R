library("rjags")
library("PhenologyBayesModeling")


calTransRoots <- function(c,b,m,d){
  root1 <- (b*m+log(sqrt(3)+2))/b
  root2 <- (b*m+log(2-sqrt(3)))/b
  return(c(root1,root2))
}
###Define functions for comparison metric calculations

##Calculate normal R^2
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

n <- 10000
#####Create data objects to be used in lots of analyses
startDay <- 1
endDay <- 365

siteData <- read.csv("GOES_Paper_Sites_FINAL.csv",header=TRUE)
#pdf("TransFigure.pdf",width=60,height=30)
#par(mfrow=c(2,3))
iseq <- c(seq(1,nrow(siteData)),seq(1,7),seq(9,15))
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


  inputFileName <- paste(siteName,"_",yr,"_GOES_NDVI_overallFilter_varBurn.RData",sep="")
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
  inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_NDVI_EIV_varBurn.RData",sep="")
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
  
  inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_MODIS_DQF_EVI_EIV_varBurn.RData",sep="")
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
  
  # #PC:
  inputFileName <- paste(siteName,"_",startDate,"_",endDate,"_PC_EIV_varBurn_SF.RData",sep="")
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
##Correct to 365
# GOES.TranFs1[GOES.TranFs1>365] <- 365
# GOES.TranFs[GOES.TranFs>365] <- 365
# MODIS.N.TranFs1[MODIS.N.TranFs1>365] <- 365
# MODIS.E.TranFs1[MODIS.E.TranFs1>365] <- 365
# MODIS.N.TranFs[MODIS.N.TranFs>365] <- 365
# MODIS.E.TranFs[MODIS.E.TranFs>365] <- 365
# PC.TranFs1[PC.TranFs1>365] <- 365
# PC.TranFs[PC.TranFs>365] <- 365

##Switch autumn TranFs and TranFs2
# GOES.TranFs1Temp <- GOES.TranFs1
# PC.TranFs1Temp <- PC.TranFs1
# MODIS.N.TranFs1Temp <- MODIS.N.TranFs1
# MODIS.E.TranFs1Temp <- MODIS.E.TranFs1
# 
# GOES.TranFs1 <- GOES.TranFs2
# PC.TranFs1 <- PC.TranFs2
# MODIS.N.TranFs1 <- MODIS.N.TranFs2
# MODIS.E.TranFs1 <- MODIS.E.TranFs2
# 
# GOES.TranFs2 <- GOES.TranFs1Temp
# PC.TranFs2 <- PC.TranFs1Temp
# MODIS.N.TranFs2 <- MODIS.N.TranFs1Temp
# MODIS.E.TranFs2 <- MODIS.E.TranFs1Temp


###Calculate Bayesian R^2
##Spring PC and GOES

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










##Identify for each site what the bias was for each transition date #First time a 0 appears in both columns

# 
# biasFiles <- dir(pattern="PhenologyBias")
# outputData <- matrix(nrow=length(iseq),ncol=length(biasFiles))
# 
# colNamesList <- character()
# 
# 
# 
# for(f in 1:length(biasFiles)){
#   newName <- paste(strsplit(biasFiles[f],"_")[[1]][2],"_",strsplit(strsplit(biasFiles[f],"_")[[1]][4],"[.]")[[1]][1],sep="")
#   print(newName)
#   colNamesList <- c(colNamesList,newName)
#   
#   inData <- read.csv(biasFiles[f],header=FALSE)
#   lessCol <- seq(2,ncol(inData),2)
#   greatCol <- seq(3,ncol(inData),2)
#   for(i in 2:(length(iseq)+1)){
#     lessBias <- min(which(inData[i,lessCol]==0))
#     greatBias <- min(which(inData[i,greatCol]==0))
#     if(lessBias==Inf){
#       lessBias <- NA
#     }
#     if(greatBias==Inf){
#       greatBias <- NA
#     }
#     if(!is.na(greatBias)&&!is.na(lessBias)){
#       
#       if(lessBias>greatBias){
#         bias <- "-" ##The bias is "<" meaning that it is < PC, earlier (-)
#       }else if(greatBias>lessBias){
#         bias <- "+"
#       }
#       else{
#         bias <- "0"
#       }
#       maxBias <- max(lessBias,greatBias)
#       if(bias!="0"){
#         bias <- paste(bias,(maxBias-1),sep="")
#       }
#     }
#     else{
#       bias <- NA
#     }
#     outputData[(i-1),f] <- as.numeric(bias)
#     #print(bias)
#     #print(maxBias)
#   }
#   
# }
# colnames(outputData) <- colNamesList
# 
# outputData <- rbind(outputData,colMeans(outputData,na.rm=TRUE))
# rownames(outputData) <- c(as.character(siteData$siteName[iseq]),"Average")
# write.csv(file="PhenologyBias_Summary.csv",outputData,col.names = TRUE,row.names = TRUE)
# 
# ####Root Mean Square Error (RMSE)
# calRMSE <- function(obs,preds){
#   sums <- 0
#   for(i in 1:length(obs)){
#     sums <- sums+(preds[i]-obs[i])**2
#   }
#   return(sqrt(sums/length(obs)))
# }
# 
# 
# 
# 
# 
# xMatrix <- PC.TranSs
# yMatrix <- MODIS.N.TranSs
# 
# R2.vals <- numeric()
# for(i in 1:10000){
#   obys <- yMatrix[,i]
#   prys <- xMatrix[,i]
#   R2.vals <- c(R2.vals,calR2(obys=obys,prys=prys))
#   
# }
# 
# hist(R2.vals)
# range(R2.vals)
# mean(R2.vals)
# quantile(R2.vals,c(0.025,0.5,0.975))
# 
# #mean(mu1 - mu2 < 0)
# bias1 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias1 <- c(bias1,mean(xMatrix[i,]-yMatrix>daylimit))
# }
# #bias1
# 
# bias2 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias2 <- c(bias2,mean(-xMatrix[i,]+yMatrix>daylimit))
# }
# #bias2
# spring.PC.MODIS.N <- cbind(bias1,bias2)
# spring.PC.MODIS.N <- cbind(spring.PC.MODIS.N,spring.PC.MODIS.N[,1]>0.95,spring.PC.MODIS.N[,2]>0.95)
# 
# spring.PC.GOES
# spring.PC.MODIS.N
# 
# #################################
# daylimit <- 1
# xMatrix <- PC.TranSs
# yMatrix <- GOES.TranSs
# 
# #mean(mu1 - mu2 < 0)
# bias1 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias1 <- c(bias1,mean(xMatrix[i,]-yMatrix[i,]>daylimit))
# }
# #bias1
# 
# bias2 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias2 <- c(bias2,mean(-xMatrix[i,]+yMatrix[i,]>daylimit))
# }
# #bias2
# spring.PC.GOES <- cbind(bias1,bias2)
# spring.PC.GOES <- cbind(spring.PC.GOES,spring.PC.GOES[,1]>0.95,spring.PC.GOES[,2]>0.95)
# 
# 
# xMatrix <- PC.TranSs
# yMatrix <- MODIS.N.TranSs
# 
# #mean(mu1 - mu2 < 0)
# bias1 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias1 <- c(bias1,mean(xMatrix[i,]-yMatrix[i,]>daylimit))
# }
# #bias1
# 
# bias2 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias2 <- c(bias2,mean(-xMatrix[i,]+yMatrix[i,]>daylimit))
# }
# #bias2
# spring.PC.MODIS.N <- cbind(bias1,bias2)
# spring.PC.MODIS.N <- cbind(spring.PC.MODIS.N,spring.PC.MODIS.N[,1]>0.95,spring.PC.MODIS.N[,2]>0.95)
# 
# spring.PC.GOES
# spring.PC.MODIS.N
# 
# ########
# daylimit <- 4
# xMatrix <- PC.TranSs1
# yMatrix <- GOES.TranSs1
# 
# #mean(mu1 - mu2 < 0)
# bias1 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias1 <- c(bias1,mean(xMatrix[i,]-yMatrix[i,]>daylimit))
# }
# #bias1
# 
# bias2 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias2 <- c(bias2,mean(-xMatrix[i,]+yMatrix[i,]>daylimit))
# }
# #bias2
# spring.PC.GOES <- cbind(spring.PC.GOES,bias1,bias2)
# #spring.PC.GOES <- cbind(spring.PC.GOES,spring.PC.GOES[,1]>0.95,spring.PC.GOES[,2]>0.95)
# 
# 
# xMatrix <- PC.TranSs1
# yMatrix <- MODIS.N.TranSs1
# 
# #mean(mu1 - mu2 < 0)
# bias1 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias1 <- c(bias1,mean(xMatrix[i,]-yMatrix[i,]>daylimit))
# }
# #bias1
# 
# bias2 <- numeric()
# for(i in 1:nrow(xMatrix)){
#   bias2 <- c(bias2,mean(-xMatrix[i,]+yMatrix[i,]>daylimit))
# }
# #bias2
# spring.PC.MODIS.N <- cbind(spring.PC.MODIS.N,bias1,bias2)
# #spring.PC.MODIS.N <- cbind(spring.PC.MODIS.N,spring.PC.MODIS.N[,1]>0.95,spring.PC.MODIS.N[,2]>0.95)
# 
# 
# ##Higher
# colnames(spring.PC.GOES) <- c("PC 1 Day","GOES 1 Day","PC 2 Day","GOES 2 Day","PC 3 Day","GOES 3 Day","PC 4 Day","GOES 4 Day")
# colnames(spring.PC.MODIS.N) <- c("PC 1 Day","MODIS.N 1 Day","PC 2 Day","MODIS.N 2 Day","PC 3 Day","MODIS.N 3 Day","PC 4 Day","MODIS.N 4 Day")
# 
# rbind(spring.PC.GOES,colSums(spring.PC.GOES))
# 
# rbind(spring.PC.MODIS.N,colSums(spring.PC.MODIS.N))
# 
# plot(rowMeans(xMatrix),rowMeans(yMatrix),pch=20)
# abline(0,1,col="red")
# 
# 
# 
# 



