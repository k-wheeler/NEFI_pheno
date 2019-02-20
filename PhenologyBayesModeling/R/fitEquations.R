##' Compute the logistic curve for phenology values from Zhang et al. (2003)
##'
##' @param Tran  DOY of transition date
##' @param b b parameter (from Zhang et al., 2003)
##' @param c c parameter (from Zhang et al., 2003)
##' @param d d parameter (from Zhang et al., 2003)
##' @param xseq the sequence of x values to compute the logistic curve for
##' @export
pheno.logistic <- function(Tran,b,c,d,xseq){
  return(c/(1+exp(b*(xseq-Tran)))+d)
}

##' Rescales the values to be between 0 and 1
##'
##' @param c c parameter (from Zhang et al., 2003)
##' @param d d parameter (from Zhang et al., 2003)
##' @param yseq the sequence of y values to compute the logistic curve for
##' @export
rescale <- function(c,d,yseq){
  return((yseq-c-d)/c+1)
}

##' Compute the logistic curve for phenology values from Zhang et al. (2003)
##'
##' @param Tran  DOY of transition date
##' @param b b parameter from shrubland process model
##' @param c c parameter from shrubland process model
##' @param d d parameter from shrubland process model
##' @param k k parameter from shrubland process model
##' @param r r parameter from shrubland process model
##' @param xseq the sequence of x values to compute the shrubland curve for
##' @export
shrublandYvals <- function(Tran,b,c,d,k,r,xseq){
  bk <- which(xseq==round(k,digits=0))

  if(k<110){
    bk <- 1
  }
  greenup <- pheno.logistic(Tran=Tran,b=b,c=c,d=d,xseq[1:bk])
  gd.xseq <- xseq[(bk+1):length(xseq)]-k
  greendown <- c*exp(r*gd.xseq)+d
  return(c(greenup,greendown))
}

##' Compute the double logistic curve (autumn then spring) for phenology values
##'
##' @param TranF  DOY of autumn/fall transition date
##' @param TranS  DOY of spring transition date
##' @param bF autumn/fall b parameter (from Zhang et al., 2003)
##' @param bS spring b parameter (from Zhang et al., 2003)
##' @param c c parameter (from Zhang et al., 2003)
##' @param d d parameter (from Zhang et al., 2003)
##' @param k k (changepoint) parameter
##' @param xseq the sequence of x values to compute the double logistic curve for
##' @param seasonOrder the order of seasons FS or SF
##' @export
deciduousYvals <- function(TranF,bF,TranS,bS,c,d,k,xseq,seasonOrder="FS"){
  bk <- which(xseq==round(k,digits=0))
  if(seasonOrder=="FS"){
    greendown <- pheno.logistic(Tran=TranF,b=bF,c=c,d=d,xseq[1:bk])
    gu.xseq <- xseq[(bk+1):length(xseq)]
    greenup <- pheno.logistic(Tran=TranS,b=bS,c=c,d=d,gu.xseq)
    return(c(greendown,greenup))
  }
  else if(seasonOrder=="SF"){
    greenup <- pheno.logistic(Tran=TranS,b=bS,c=c,d=d,xseq[1:bk])
    gd.xseq <- xseq[(bk+1):length(xseq)]
    greendown <- pheno.logistic(Tran=TranF,b=bF,c=c,d=d,gd.xseq)
    return(c(greenup,greendown))
  }

}

##' Compute the gaussian curve values
##'
##' @param a a parameter
##' @param b b parameter
##' @param c c parameter
##' @param d d parameter
##' @param xseq the sequence of x values to compute the double gaussian curve for
##' @export
gaussianCurve <- function(a,b,c,d,xseq){
  return(a*exp((-1*(xseq-b)**2)/(2*c**2))+d)
}
