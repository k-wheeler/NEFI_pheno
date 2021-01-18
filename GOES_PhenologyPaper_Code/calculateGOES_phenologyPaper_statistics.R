####Root Mean Square Error (RMSE)
calRMSE <- function(obs,preds){
  sums <- 0
  for(i in 1:length(obs)){
    sums <- sums+(preds[i]-obs[i])**2
  }
  return(sqrt(sums/length(obs)))
}


##Will need to loop over the combinations (maybe start with one basic: GOES vs PC spring)
# IDs <- c("PC_GOES_Spring50","PC_MODISN_Spring50","PC_MODISE_Spring50","MODISN_GOES_Spring50","MODISE_GOES_Spring50","MODISE_MODISN_Spring50",
#          "PC_GOES_Fall50","PC_MODISN_Fall50","PC_MODISE_Fall50","MODISN_GOES_Fall50","MODISE_GOES_Fall50","MODISE_MODISN_Fall50",
#          "PC_GOES_SpringTran1","PC_MODISN_SpringTran1","PC_MODISE_SpringTran1","MODISN_GOES_SpringTran1","MODISE_GOES_SpringTran1","MODISE_MODISN_SpringTran1",
#          "PC_GOES_FallTran1","PC_MODISN_FallTran1","PC_MODISE_FallTran1","MODISN_GOES_FallTran1","MODISN_GOES_FallTran1","MODISE_MODISN_FallTran1",
#          "PC_GOES_SpringTran2","PC_MODISN_SpringTran2","PC_MODISE_SpringTran2","MODISN_GOES_SpringTran2","MODISE_GOES_SpringTran2","MODISE_MODISN_SpringTran2",
#          "PC_GOES_FallTran2","PC_MODISN_FallTran2","PC_MODISE_FallTran2","MODISN_GOES_FallTran2","MODISE_GOES_FallTran2","MODISE_MODISN_FallTran2")

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


##Load and add bias data to output table
biasFiles <- dir(pattern="PhenologyBias")
biases <- character()
for(i in 1:length(IDs)){
  fileName <- paste("PhenologyBias_",IDs[i],".csv",sep="")
  dat <- read.csv(fileName,header=TRUE)
  
  #newBias <- paste(round(dat[30,3],digits=2)," (",round(dat[30,2],digits=2),", ",round(dat[30,4],digits=2),")",sep="")
  newBias <- paste(round(dat[30,3],digits=2)," (",round((dat[30,4]-dat[30,2])/2,digits=2),")",sep="")
  biases <- c(biases,newBias)
}
output <- cbind(IDs,as.numeric(round(R2s,digits=2)),as.numeric(round(RMSEs,digits=2)),biases)
colnames(output) <- c("ID","R^2","RMSE","Average Bias (95% CI)")

write.table(output,file="GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_FitStatisticsUpdated.csv",row.names = FALSE,col.names = colnames(output),sep=",")

##p-values: calculate the predicted quantile of each observed data point (aka the p-value); 
#histogram of Bayesian p-values should have a uniform distribution

##From median estimate of PC transition dates

# PC.TranFs1med <- apply(PC.TranFs1,1,quantile,0.5)
# PC.TranFs2med <- apply(PC.TranFs2,1,quantile,0.5)
# PC.TranSs1med <- apply(PC.TranSs1,1,quantile,0.5)
# PC.TranSs2med <- apply(PC.TranSs2,1,quantile,0.5)

# plotPValues <- function(fileStr,sourceName,season,TranFs1,TranFs2,TranSs1,TranSs2,TranSs,TranFs){
#   p.valsF50 <- matrix(nrow=0,ncol=3)
#   p.valsS50 <- matrix(nrow=0,ncol=3)
#   p.valsF1 <- matrix(nrow=0,ncol=3)
#   p.valsF2 <- matrix(nrow=0,ncol=3)
#   p.valsS1 <- matrix(nrow=0,ncol=3)
#   p.valsS2 <- matrix(nrow=0,ncol=3)
#   
#   for(i in 1:length(iseq)){
#     siteName <- as.character(siteData$siteName[iseq[i]])
#     print(siteName)
#     inputFileName <- paste(siteName,fileStr,sep="")
#     print(inputFileName)
#     #load(inputFileName)
#     
#     #var.mat <- as.matrix(var.burn)
#     
#     #print(colnames(var.mat))
#     #GOES.TranF<-var.mat[,1]
#     #GOES.TranS<-var.mat[,2]
#     
#     GOES.TranF1 <- TranFs1[i,]
#     GOES.TranF2 <- TranFs2[i,]
#     GOES.TranS1 <- TranSs1[i,]
#     GOES.TranS2 <- TranSs2[i,]
#     GOES.TranF <- TranFs[i,]
#     GOES.TranS <- TranSs[i,]
#     
#     p.valsF50new <- numeric()
#     p.valsS50new <- numeric()
#     p.valsF1new <- numeric()
#     p.valsF2new <- numeric()
#     p.valsS1new <- numeric()
#     p.valsS2new <- numeric()
#     
#     for(g in 1:10000){
#       
#       p.valsF50new <- c(p.valsF50new,sum(PC.TranFs[i,g] > GOES.TranF)/nrow(var.mat))
#       p.valsS50new <- c(p.valsS50new,sum(PC.TranSs[i,g] > GOES.TranS)/nrow(var.mat))
#       
#       p.valsF1new <- c(p.valsF1new,sum(PC.TranFs1[i,g] > GOES.TranF1)/length(GOES.TranF1))
#       p.valsS1new <- c(p.valsS1new,sum(PC.TranSs1[i,g] > GOES.TranS1)/length(GOES.TranS1))
#       
#       p.valsF2new <- c(p.valsF2new,sum(PC.TranFs2[i,g] > GOES.TranF2)/length(GOES.TranF2))
#       p.valsS2new <- c(p.valsS2new,sum(PC.TranSs2[i,g] > GOES.TranS2)/length(GOES.TranS2))
#       
#       # 
#       #     p.valsF50 <- c(p.valsF50,sum(quantile(PC.TranFs[i,],0.50) > GOES.TranF)/nrow(var.mat))
#       #     p.valsS50 <- c(p.valsS50,sum(quantile(PC.TranSs[i,],0.50) > GOES.TranS)/nrow(var.mat))
#       #     
#       #     p.valsF1 <- c(p.valsF1,sum(quantile(PC.TranFs1med[i],0.50) > GOES.TranF1)/length(GOES.TranF1))
#       #     p.valsS1 <- c(p.valsS1,sum(quantile(PC.TranSs1med[i],0.50) > GOES.TranS1)/length(GOES.TranS1))
#       #     
#       #     p.valsF2 <- c(p.valsF2,sum(quantile(PC.TranFs2med[i],0.50) > GOES.TranF2)/length(GOES.TranF2))
#       #     p.valsS2 <- c(p.valsS2,sum(quantile(PC.TranSs2med[i],0.50) > GOES.TranS2)/length(GOES.TranS2))
#     }
#     p.valsF50 <- rbind(p.valsF50,quantile(p.valsF50new,c(0.025,0.5,0.975)))
#     p.valsS50 <- rbind(p.valsS50,quantile(p.valsS50new,c(0.025,0.5,0.975)))
#     
#     p.valsF1 <- rbind(p.valsF1,quantile(p.valsF1new,c(0.025,0.5,0.975)))
#     p.valsS1 <- rbind(p.valsS1,quantile(p.valsS1new,c(0.025,0.5,0.975)))
#     
#     p.valsF2 <- rbind(p.valsF2,quantile(p.valsF2new,c(0.025,0.5,0.975)))
#     p.valsS2 <- rbind(p.valsS2,quantile(p.valsS2new,c(0.025,0.5,0.975)))
#     #print(quantile(p.valsF50,c(0.025,0.5,0.975)))
#     #print(quantile(p.valsS50,c(0.025,0.5,0.975)))
#     
#   }
#   return(list(p.valsF50=p.valsF50,p.valsF1=p.valsF1,p.valsF2=p.valsF2,p.valsS50=p.valsS50,p.valsS1=p.valsS1,p.valsS2=p.valsS2))
#   #sourceName <- "MODIS EVI"
#   if(season=="fall"){
#     #h <- hist(p.valsF1, breaks = 10, plot=FALSE)
#     #h$counts=h$counts/sum(h$counts)
#     plot(density(p.valsF1),xlab="P-values",main=paste(sourceName, "Autumn Transition 1"),xlim=c(0,1)) 
#     
#     #h <- hist(p.valsF50, breaks = 10, plot=FALSE)
#     #h$counts=h$counts/sum(h$counts)
#     plot(density(p.valsF50),xlab="P-values",main=paste(sourceName, "Fall 50%"),xlim=c(0,1)) 
#     
#     #h <- hist(p.valsF2, breaks = 10, plot=FALSE)
#     #h$counts=h$counts/sum(h$counts)
#     plot(density(p.valsF2),xlab="P-values",main=paste(sourceName, "Autumn Transition 2"),xlim=c(0,1)) 
#     
#   }else{
#     #h <- hist(p.valsS1, breaks = 10, plot=FALSE)
#     #h$counts=h$counts/sum(h$counts)
#     plot(density(p.valsS1),xlab="P-values",main=paste(sourceName, "Spring Transition 1"),xlim=c(0,1)) 
#     
#     #h <- hist(p.valsS50, breaks = 10, plot=FALSE)
#     #h$counts=h$counts/sum(h$counts)
#     plot(density(p.valsS50),xlab="P-values",main=paste(sourceName, "Spring 50%"),xlim=c(0,1)) 
#     
#     #h <- hist(p.valsS2, breaks = 10, plot=FALSE)
#     #h$counts=h$counts/sum(h$counts)
#     plot(density(p.valsS2),xlab="P-values",main=paste(sourceName, "Spring Transition 2"),xlim=c(0,1)) 
#   }
# }
# GOES.pDat <- plotPValues(fileStr = "_overall6_varBurn.RData",season="spring",sourceName="GOES",TranFs1=GOES.TranFs1,TranFs2=GOES.TranFs2,TranSs1=GOES.TranSs1,TranSs2=GOES.TranSs2,TranFs=GOES.TranFs,TranSs=GOES.TranSs)
# save(file="GOES.pDat.R",GOES.pDat)
# MODISN.pDat <- plotPValues(fileStr = "_182_546_MODIS_DQF_NDVI_varBurn.RData",season="spring",sourceName="MODIS NDVI",TranFs1=MODIS.N.TranFs1,TranFs2=MODIS.N.TranFs2,TranSs1=MODIS.N.TranSs1,TranSs2=MODIS.N.TranSs2,TranFs=MODIS.N.TranFs,TranSs=MODIS.N.TranSs)             
# save(file="MODISN.pDat.R",MODISN.pDat)
# MODISE.pDat <- plotPValues(fileStr = "_182_546_MODIS_DQF_EVI_varBurn.RData",season="spring",sourceName="MODIS EVI",TranFs1=MODIS.E.TranFs1,TranFs2=MODIS.E.TranFs2,TranSs1=MODIS.E.TranSs1,TranSs2=MODIS.E.TranSs2,TranFs=MODIS.E.TranFs,TranSs=MODIS.E.TranSs)
# save(file="MODISE.pDat.R",MODISE.pDat)
# load("MODISN.pDat.R")
# 
# #pdf("TransComparison_pValues.pdf",width=15,height=10)
# jpeg("GOES_PaperWriting/GOES_PhenologyPaper/GOES_Phenology_Paper_Transition_pVales_RAW.jpeg",width=6,height=3,units = "in",res=1000)  
# 
# par(mfrow=c(2,3),mai=c(0.3,0.1,0.1,0.1))
# plot(density(MODISE.pDat$p.valsS1[,2]),col=MODIS.E.col,xlim=c(0,1),main="",lwd=2,xlab="",frame.plot=FALSE,yaxt="n",ylab="",cex.lab=2)
# #ciEnvelope(density(GOES.pDat$p.valsS1[,1])$x,density(GOES.pDat$p.valsS1[,1])$y,density(GOES.pDat$p.valsS1[,3])$y,col=adjustcolor(GOES.col,0.8))
# #ciEnvelope(density(MODISN.pDat$p.valsS1[,1])$x,density(MODISN.pDat$p.valsS1[,1])$y,density(MODISN.pDat$p.valsS1[,3])$y,col=adjustcolor(MODIS.N.col,0.6))
# 
# #ciEnvelope(density(MODISE.pDat$p.valsS1[,1])$x,density(MODISE.pDat$p.valsS1[,1])$y,density(MODISE.pDat$p.valsS1[,3])$y,col=adjustcolor(MODIS.E.col,0.4))
# lines(density(MODISE.pDat$p.valsS1[,2]),col=MODIS.E.col,lwd=2)
# lines(density(MODISN.pDat$p.valsS1[,2]),col=MODIS.N.col,lwd=2)
# lines(density(GOES.pDat$p.valsS1[,2]),col=GOES.col,lwd=2)
# legend("topleft",c("GOES","MODIS NDVI","MODIS EVI"),lwd=c(2,2,2),col=c(GOES.col,MODIS.N.col,MODIS.E.col))
# 
# 
# plot(density(GOES.pDat$p.valsS50[,2]),col=GOES.col,xlim=c(0,1),main="",lwd=2,xlab="",frame.plot=FALSE,yaxt="n",ylab="",cex.lab=2)
# lines(density(MODISE.pDat$p.valsS50[,2]),col=MODIS.E.col,lwd=2)
# lines(density(MODISN.pDat$p.valsS50[,2]),col=MODIS.N.col,lwd=2)
# 
# plot(density(GOES.pDat$p.valsS2[,2]),col=GOES.col,xlim=c(0,1),main="",lwd=2,xlab="",frame.plot=FALSE,yaxt="n",ylab="")
# lines(density(MODISE.pDat$p.valsS2[,2]),col=MODIS.E.col,lwd=2)
# lines(density(MODISN.pDat$p.valsS2[,2]),col=MODIS.N.col,lwd=2)
# 
# plot(density(MODISE.pDat$p.valsF2[,2]),col=MODIS.E.col,xlim=c(0,1),main="",lwd=2,xlab="",frame.plot=FALSE,yaxt="n",ylab="",cex.lab=2) #Fall trans1 and trans2 are different
# lines(density(GOES.pDat$p.valsF2[,2]),col=GOES.col,lwd=2)
# lines(density(MODISN.pDat$p.valsF2[,2]),col=MODIS.N.col,lwd=2)
# 
# plot(density(MODISE.pDat$p.valsF50[,2]),col=MODIS.E.col,xlim=c(0,1),main="",lwd=2,xlab="",frame.plot=FALSE,yaxt="n",ylab="",cex.lab=2)
# lines(density(GOES.pDat$p.valsF50[,2]),col=GOES.col,lwd=2)
# lines(density(MODISN.pDat$p.valsF50[,2]),col=MODIS.N.col,lwd=2)
# 
# plot(density(GOES.pDat$p.valsF1[,2]),col=GOES.col,xlim=c(0,1),main="",lwd=2,xlab="",frame.plot=FALSE,yaxt="n",ylab="",cex.lab=2)
# lines(density(MODISN.pDat$p.valsF1[,2]),col=MODIS.N.col,lwd=2)
# lines(density(MODISE.pDat$p.valsF1[,2]),col=MODIS.E.col,lwd=2)
# 
# dev.off()
# ############
# 
# 
# 
# 





