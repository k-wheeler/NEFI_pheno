##' Creates an output summary (not finalized)
##'
##' @param siteFileName The desired output file name
##' @param siteName iseq The row numbers of the sites in the data file to include in the summary
##' @param startDay The start day
##' @param endDay The end day
createOutputSummary <- function(siteFileName,iseq,startDay,endDay){
  siteData <- read.csv(siteFileName,header=FALSE)
  xseq <- seq(startDay,endDay,1)
  for(i in iseq){
    siteName <- as.character(siteData[i,1])
    print(siteName)
    URL <- as.character(siteData[i,4])
    PFT <- as.character(siteData[i,5])

    #GOES
    data.GOES = GOES_data(siteName=siteName,startDay = startDay,endDay = endDay)
    inputFileName <- paste(siteName,"_GOES_varBurnINITS.RData",sep="")
    load(inputFileName)
    print("GOES Summary:")
    print(summary(GOES.md.out))
    # GOES.c <- var.mat[,3]
    # GOES.d <- var.mat[,4]
    # if(PFT=="DB"){
    #   GOES.a <- var.mat[,1]
    #   GOES.b <- var.mat[,2]
    #   GOES.trans <- mean(-(GOES.a/GOES.b))
    # }
    # else if(PFT=="SH"){
    #   GOES.trans <- mean(var.mat[,5])
    # }
    # ci.GOES <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)

    #MODIS:
    data.MODIS = MODIS_data(siteName=siteName)
    inputFileName <- paste(siteName,"_MODIS_varBurnINITS.RData",sep="")
    load(inputFileName)
    print("MODIS Summary")
    print(summary(MODIS.md.out))

    # var.mat<-as.matrix(MODIS.md.out)
    # MODIS.c <- var.mat[,3]
    # MODIS.d <- var.mat[,4]
    #
    # if(PFT=="DB"){
    #   MODIS.a <- var.mat[,1]
    #   MODIS.b <- var.mat[,2]
    #   MODIS.trans <- mean(-(MODIS.a/MODIS.b))
    # }
    # else if(PFT=="SH"){
    #   MODIS.trans <- mean(var.mat[,5])
    # }
    # ci.MODIS <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)
    # print("MODIS Done")

    #PC:
    data.PC = PC_data(URL=URL,startDay=startDay,endDay=endDay)
    inputFileName <- paste(siteName,"_PC_varBurnINITS.RData",sep="")
    load(inputFileName)
    print("PC Summary:")
    print(summary(PC.md.out))

    # var.mat<-as.matrix(PC.md.out)
    # PC.c <- var.mat[,3]
    # PC.d <- var.mat[,4]
    # if(PFT=="DB"){
    #   PC.a <- var.mat[,1]
    #   PC.b <- var.mat[,2]
    #   PC.trans <- mean(-(PC.a/PC.b))
    # }
    # else if(PFT=="SH"){
    #   PC.trans <- mean(var.mat[,5])
    # }
    # ci.PC <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)
    # print("PC Done")
  }
}
