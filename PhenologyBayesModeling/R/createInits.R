##' Function to create initial conditions
##'
##' @param data The data in the form of a list with items x and y
##' @param PFT Plant functional type ("SH" or "DB")
##' @export
createInits <- function(data,PFT){
  x <- data$x
  y <- data$y
  d <- mean(sort(y)[1:10],na.rm=TRUE)
  c <- mean(sort(y,decreasing = TRUE)[1:5],na.rm=TRUE)-d
  init.mus <- list(d=d,c=c)
  if(PFT=="SH"){
    init.mus$k <- x[which(y==max(y,na.rm = TRUE))]
  }
  return(init.mus)
}


