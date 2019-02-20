##' Create the credible interval from MCMC output matrix
##'
##' @param PFT The plant functional type (DB or SH)
##' @param var.mat The MCMC output as a matrix
##' @param xseq The sequence of x values you want to calculate the CI for
##' @param doRescale Boolean for if you want to rescale the CI to be between 0 and 1 (default=True)
##' @param seasonOrder the order of seasons FS or SF
##' @export
createCI <- function(PFT,var.mat,xseq,doRescale=TRUE,seasonOrder="FS"){
  ycred <- matrix(0,nrow=10000,ncol=length(xseq))
  if(PFT=="DB"){
    TranF<-var.mat[,1]
    rndNums <- sample(1:length(TranF),10000,replace=T)
    TranF <- TranF[rndNums]
    TranS<-var.mat[rndNums,2]
    bF <- var.mat[rndNums,3]
    bS <- var.mat[rndNums,4]
    c <- var.mat[rndNums,5]
    d <- var.mat[rndNums,6]
    k <- var.mat[rndNums,7]
    for(g in 1:10000){
      if(doRescale){
        Ey <- rescale(c=c[g],d=d[g],yseq=deciduousYvals(TranS=TranS[g],bS=bS[g],TranF=TranF[g],bF=bF[g],c=c[g],d=d[g],k=k[g],xseq=xseq,seasonOrder = seasonOrder))
      }
      else{
      Ey <- deciduousYvals(TranS=TranS[g],bS=bS[g],TranF=TranF[g],bF=bF[g],c=c[g],d=d[g],k=k[g],xseq=xseq,seasonOrder=seasonOrder)
      }
      ycred[g,] <- Ey
    }
  }
  else if(PFT=="SH"){
    Tran<-var.mat[,1]
    rndNums <- sample(1:length(Tran),10000,replace=T)
    Tran <- Tran[rndNums]
    b<-var.mat[rndNums,2]
    c <- var.mat[rndNums,3]
    d <- var.mat[rndNums,4]
    k <- var.mat[rndNums,5]
    r <- var.mat[rndNums,7]

    for(g in 1:10000){
      if(doRescale){
      Ey <- rescale(c=c[g],d=d[g],yseq=shrublandYvals(Tran=Tran[g],b=b[g],c=c[g],d=d[g],k=k[g],r=r[g],xseq=xseq))
      }
      else{
        Ey <- shrublandYvals(Tran=Tran[g],b=b[g],c=c[g],d=d[g],k=k[g],r=r[g],xseq=xseq)
      }
      ycred[g,] <- Ey
    }
  }
  return(apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE))
}


