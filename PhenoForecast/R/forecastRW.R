#' Random walk forecast step
#'
#' @param IC Initial conditions
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
forecastRW <- function(IC,Q=0,n,NT){
  #for(i in 2:n){
  #  x[i]~dnorm(x[i-1],p.proc)
  #}
  x <- matrix(NA,n,NT)
  Xprev <- IC
  for(t in 1:NT){
    x[,t] <- rnorm(n,Xprev,Q)
    Xprev <- x[,t]
  }
  return(x)
}
