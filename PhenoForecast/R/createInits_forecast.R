#' Create inits based on previous out.burn
#'
#' @param out.burn
#' @param variableNames
#'
#' @return
#' @export
#'
createInits_forecast <- function(out.burn,variableNames){
  inits <- list()
  out.mat <- data.frame(as.matrix(out.burn))
  inits$p.proc <- mean(out.mat$p.proc)
  inits$fallLength <- mean(out.mat$fallLength)
  inits$MOF <- mean(out.mat$MOF)
  inits$sSlope <- mean(out.mat$sSlope)
  if("p.PC" %in% variableNames){
    inits$p.PC <- mean(out.mat$p.PC)
  }
  if("p.MN" %in% variableNames){
    inits$p.MN <- mean(out.mat$p.MN)
  }
  if("p.ME" %in% variableNames){
    inits$p.ME <- mean(out.mat$p.ME)
  }
  if("baseTemp" %in% variableNames){
    inits$baseTemp <- mean(out.mat$baseTemp)
  }
  return(inits)
}
