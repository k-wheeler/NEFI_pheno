#' Logistic with covariates forecast step
#'
#' @param IC Initial conditions
#' @param trans The parameter trans
#' @param b1 The parameter b1
#' @param Sf The covariate Sf (cummulative average temperatures)
#' @param Q Process error (default = 0 for deterministic runs)
#' @param n Size of Monte Carlo ensemble
#' @param NT number of days for forecast
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#'
#' @return
#' @export
#'
#' @examples
forecastLogCov <- function(IC,trans,b1,Sf,Q=0,n,NT){
  b0 <-  b1 * -1 * trans
  x <- matrix(NA,n,NT)
  Xprev <- IC
  ##For the first round:
  rl <- b1 * Sf[,1] + b0
  print(rl)
  if(rl<0){rl=0}
  r <- matrix(NA,n,NT)
  rPrev <- rl
  t <- 1 #Done for the first time step
  mu <- Xprev+ rl * Xprev * (1-Xprev)  ## latent process
  xl <- rnorm(n,mu,Q)
  xNew <- numeric()
  for(i in 1:length(xl)){
    xNew <- c(xNew,max(0, min(1,xl[i])))
  }
  x[,t] <- xNew ## trunate normal process error
  r[,t] <- rl
  Xprev <- x[,t]
  rPrev <- r[,t]

  for(t in 2:NT){
    rl <- b1 * Sf[,t] + b0
    rNew <- numeric()
    for(i in 1:length(rl)){
      if(rPrev[i]<0){
        rPrev[i] <- 0
      }
      if(rl[i]<rPrev[i]){
        rNew <- c(rNew,rPrev[i])
      }else{
        rNew <- c(rNew,rl[i])
      }
    }
    mu <- Xprev+ rNew * Xprev * (1-Xprev)  ## latent process
    xl <- rnorm(n,mu,Q)
    xNew <- numeric()
    for(i in 1:length(xl)){
      xNew <- c(xNew,max(0, min(1,xl[i])))
    }
    x[,t] <- xNew ## trunate normal process error
    r[,t] <- rNew
    Xprev <- x[,t]
    rPrev <- r[,t]
  }
  return(x)
}
