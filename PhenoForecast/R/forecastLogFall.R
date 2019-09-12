#' Basic logistic forecast step
#'
#' @param IC Initial conditions
#' @param r The parameter r
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
forecastLogFall <- function(IC,r,Q=0,n,NT){
  x <- matrix(NA,n,NT)
  Xprev <- IC
  for(t in 1:NT){
    mu <- Xprev - r * Xprev * (1-Xprev)  ## latent process
    xl <- rnorm(n,mu,Q)
    xNew <- numeric()
    for(i in 1:length(xl)){
      xNew <- c(xNew,max(0, min(1,xl[i])))
    }
    x[,t] <- xNew ## trunate normal process error
    Xprev <- x[,t]
  }
  return(x)
}
# xseq <- seq(182,220)
# yseq <- forecastLogFall(IC=0.99,r=0.3,n=1,NT=length(xseq))
# plot(xseq,yseq,ylim=c(0,1))
