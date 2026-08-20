###########Old JAGS functional form senescence models:
##Investigating priors on forecast model calibration
else if(vars=="b5"){ #full 
  #outputFileName <- paste(siteName,"_b5_triggerPara_calibration_varBurn.R",sep="")
  outputFileName <- paste(siteName,"_b5_triggerPara_Dtrigger_calibration_varBurn.RData",sep="")
  #outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_b5_triggerPara_calibration_varBurn.R",sep="")
  variableNames <- c("p.PC","p.proc","x","b0","b1","b2","b5","Dtrigger")
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 0
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 1
  dataFinal$b5_lower <- 0
  dataFinal$b5_upper <- 1
  #dataFinal$CDDtrigger.lower <- 0
  #dataFinal$CDDtrigger.upper <- 500
  dataFinal$Dtrigger.lower <- 10
  dataFinal$Dtrigger.upper <- 20
  
  generalModel = "
  model {
  ### Data Models for complete years
  for(yr in 1:(N-1)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N-1)){
  for(i in 2:n){
  
  ##Calculate CDD
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  #CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
  
  ##offset[i,yr] <- (baseTemp - Tair[i,yr])
  ##newCDD1[i,yr] <- ifelse(offset[i,yr]>0,offset[i,yr],0)
  
  ##newCDD2[i,yr] <- CDDs[(i-1),yr] + newCDD1[i,yr]
  ##CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,newCDD2[i,yr],CDDs[(i-1),yr])
  
  
  #xmu[i,yr] <- x[(i-1),yr] + ifelse(CDDs[i,yr]<CDDtrigger,0,(b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2) + (b5 * Tair[i,yr] * D[i,yr])))
  xmu[i,yr] <- x[(i-1),yr] + ifelse(D[i,yr]<Dtrigger,0,(b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2) + (b5 * Tair[i,yr] * D[i,yr])))
  x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
  x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
  
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a,x1.b)
  #CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  #CDDtrigger ~ dnorm(CDDtrigger.mu,CDDtrigger.prec)
  #CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
  Dtrigger ~ dunif(Dtrigger.lower,Dtrigger.upper)
  
  b0 ~ dunif(b0_lower,b0_upper)
  b1 ~ dunif(b1_lower,b1_upper)
  b2 ~ dunif(b2_lower,b2_upper)
  b5 ~ dunif(b5_lower,b5_upper)
  }
  
  "
}else if(vars=="b5_exp"){ #full 
  
  outputFileName <- paste(siteName,"_b5_exp_noTrigger_calibration_varBurn.RData",sep="")
  variableNames <- c("p.PC","p.proc","x","b1","b5")
  
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  dataFinal$b5_lower <- 0
  dataFinal$b5_upper <- 1
  
  generalModel = "
  model {
  ### Data Models for complete years
  for(yr in 1:(N-1)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N-1)){
  for(i in 2:n){
  
  ##Calculate CDD
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  xmu[i,yr] <- x[(i-1),yr] + (b1 * x[(i-1),yr]) + (b5 * Tair[i,yr] * D[i,yr])
  x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
  x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a,x1.b)
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  b1 ~ dunif(b1_lower,b1_upper)
  b5 ~ dunif(b5_lower,b5_upper)
  }
  
  "
  }else if(vars=="b5_lin"){ #full 
    
    outputFileName <- paste(siteName,"_b5_lin_noTrigger_setb0_calibration_varBurn.RData",sep="")
    #variableNames <- c("p.PC","p.proc","x","b0","b5")
    variableNames <- c("p.PC","p.proc","x","b5")
    
    #dataFinal$b0_lower <- -1
    #dataFinal$b0_upper <- 0
    dataFinal$b5_lower <- 0
    dataFinal$b5_upper <- 1
    
    generalModel = "
    model {
    ### Data Models for complete years
    for(yr in 1:(N-1)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    
    #### Process Model
    for(yr in 1:(N-1)){
    for(i in 2:n){
    
    Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
    xmu[i,yr] <- x[(i-1),yr] + b0 + (b5 * Tair[i,yr] * D[i,yr])
    x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
    x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
    }
    }
    
    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a,x1.b)
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    #b0 ~ dunif(b0_lower,b0_upper)
    b0 <- 0.50
    b5 ~ dunif(b5_lower,b5_upper)
    }
    
    "
    }else if(vars=="linear"){
      outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_linear_trigger0_calibration_varBurn.RData",sep="")
      dataFinal$b0_lower <- -1
      dataFinal$b0_upper <- 0
      
      variableNames <- c("p.PC","p.proc","x","b0")
      
      generalModel = "
      model {
      ### Data Models for complete years
      for(yr in 1:(N)){
      for(i in 1:n){
      p[i,yr] ~ dnorm(x[i,yr],p.PC)
      }
      }
      
      #### Process Model
      for(yr in 1:(N)){
      for(i in 2:n){
      
      ##Calculate CDD
      Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
      #CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
      
      xmu[i,yr] <- x[(i-1),yr] + ifelse(Tair[i,yr]<baseTemp,b0,0)
      x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
      x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
      
      }
      }
      
      #### Priors
      for(yr in 1:N){ ##Initial Conditions
      x[1,yr] ~ dbeta(x1.a,x1.b)
      #CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
      }
      p.PC ~ dgamma(s1.PC,s2.PC)
      p.proc ~ dgamma(s1.proc,s2.proc)
      
      #CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
      
      b0 ~ dunif(b0_lower,b0_upper)
      #b1 ~ dunif(b1_lower,b1_upper)
      
      }
      
      "
      
    }
}else if(vars=="noCov_changepoint"){ #no covariates and a changepoint determined by GCC for start of senescence 
  outputFileName <- paste(siteName,"noCov_GCCchangepoint_calibration_includeAllBs_varBurn.RData",sep="")
  outputFileName <- paste(siteName,"noCov_GCCchangepoint90_calibration_includeAllBs_varBurn.RData",sep="")
  outputFileName <- paste(siteName,"noCov_GCCchangepoint85Trigger_calibration_includeAllBs_varBurn.RData",sep="")
  #outputFileName <- paste(siteName,"noCov_GCCchangepoint90_calibration_onlyB0_varBurn.RData",sep="")
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 0
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 1
  slps <- numeric()
  for(yr in 1:ncol(dataFinal$p)){
    slps <- c(slps,min(as.numeric(lm(dataFinal$p[1:20,yr]~seq(1,20))$coefficients[2]),0))
  }
  dataFinal$summerDeclines <- slps
  #print(dataFinal$summerDeclines)
  #print(dataFinal$p[1:10,])
  
  variableNames <- c("p.PC","p.proc","x","b0","b1","b2")#,"CDDtrigger")
  #variableNames <- c("p.PC","x","b0")#,"b1","b2")#,"CDDtrigger")
  #variableNames <- c("changepoint")
  
  ##Options for the changepoint based on GCC:
  ##1) GCC reaches a certain value; 
  ##2) GCC reaches a certain value that is x amount less than the mean summer value
  ##3) GCC reaches a certain value that is x percent less than the mean summer value
  ##4) GCC decreases by x amount; by x percent
  ##5) Average GCC decreases by x amount; by x percent
  ##I probably need to plot all of them and see which one might work the best. But maybe for right now just 
  ##chose a basic one (#5)
  generalModel = "
  model {
  for(yr in 1:(N)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  #### Process Model
  for(yr in 1:N){
  for(i in 2:n){
  #changepoint[i,yr] <- max(ifelse((x[(i-1),yr]-x[(i-2),yr]<1.05*summerDeclines[yr]),1,0),changepoint[(i-1),yr]) #1 if sensed (more negative slope); 0 if false 
  changepoint[i,yr] <- max(ifelse(x[(i-1),yr]<0.85,1,0),changepoint[(i-1),yr])
  
  xmu[i,yr] <- x[(i-1),yr] + ifelse(changepoint[i,yr]==1, (b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)), summerDeclines[yr])
  
  #xmu[i,yr] <- x[(i-1),yr] + ifelse(x[(i-1),yr]<0.90, (b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)), summerDeclines[yr])
  #xmu[i,yr] <- x[(i-1),yr] + ifelse(x[(i-1),yr]<0.90, b0, summerDeclines[yr])
  x[i,yr] ~ dnorm(xmu[i,yr],p.proc) T(0,0.999)
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a,x1.b)
  #x[2,yr] ~ dbeta(x1.a,x1.b)
  changepoint[1,yr] <- 0
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  b0 ~ dunif(b0_lower,b0_upper)
  b1 ~ dunif(b1_lower,b1_upper)
  b2 ~ dunif(b2_lower,b2_upper)
  }
  "
  
}else if(vars=="exponential"){ #Only b1
  
  #outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_exponential_triggerPara_withTrigger_calibration_varBurn.RData",sep="")
  outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_exponential_trigger0_calibration_varBurn.RData",sep="")
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  # dataFinal$CDDtrigger.lower <- 0
  # dataFinal$CDDtrigger.upper <- 500
  
  variableNames <- c("p.PC","p.proc","x","b1")
  
  generalModel = "
  model {
  ### Data Models for complete years
  for(yr in 1:(N)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N)){
  for(i in 2:n){
  
  ##Calculate CDD
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  #CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
  
  xmu[i,yr] <- x[(i-1),yr] + ifelse(Tair[i,yr]<baseTemp,(b1 * x[(i-1),yr]),0)
  x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
  x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
  
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a,x1.b)
  #CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  #CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
  
  #b0 ~ dunif(b0_lower,b0_upper)
  b1 ~ dunif(b1_lower,b1_upper)
  #b2 ~ dunif(b2_lower,b2_upper)
  #b3 ~ dunif(b3_lower,b3_upper)
  #b4 ~ dunif(b4_lower,b4_upper)
  #b5 ~ dunif(b5_lower,b5_upper)
  }
  
  "
  
}else if(vars=="nob2"){
  #outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_nob2_setTrigger_calibration_varBurn.RData",sep="")
  outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_nob2_trigger0_calibration_varBurn.RData",sep="")
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 0
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  #dataFinal$CDDtrigger <- 200
  
  variableNames <- c("p.PC","p.proc","x","b0","b1")
  
  generalModel = "
  model {
  ### Data Models for complete years
  for(yr in 1:(N)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N)){
  for(i in 2:n){
  
  ##Calculate CDD
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  #CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
  xmu[i,yr] <- x[(i-1),yr] + ifelse(Tair[i,yr]<baseTemp,(b0 + (b1 * x[(i-1),yr])),0)
  x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
  x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
  
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a,x1.b)
  #CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  #CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
  
  b0 ~ dunif(b0_lower,b0_upper)
  b1 ~ dunif(b1_lower,b1_upper)
  
  }
  
  "
  
}else if(vars=="invExponential"){ #Only b1
  #dataFinal$CDDtrigger.lower <- 0
  #dataFinal$CDDtrigger.upper <- 500
  #dataFinal$Dtrigger.lower <- 10
  #dataFinal$Dtrigger.upper <- 20
  # inits <- list()
  # 
  # for(i in 1:nchain){
  #   inits[[i]] <- list(b0=rnorm(1,-0.15,0.01))#,CDDtrigger=rnorm(1,50,5))
  # }
  #outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_invExponential_setTrigger_calibration_varBurn.RData",sep="")
  outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_invExponential_trigger0_newPriors_newConstraint_calibration_varBurn.RData",sep="")
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 0
  #dataFinal$CDDtrigger <- 200
  
  variableNames <- c("p.PC","p.proc","x","b0")#,"Dtrigger")
  
  generalModel = "
  model {
  ### Data Models for complete years
  for(yr in 1:(N)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N)){
  for(i in 2:n){
  
  ##Calculate CDD
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  #CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
  
  #deltaX[i,yr] <- ifelse(CDDs[i,yr]>CDDtrigger,b0 + (-1 * b0 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2) + (b3 * Tair[i,yr]) + (b4 * D[i,yr]) + (b5 * Tair[i,yr] * D[i,yr]),0)  
  
  #xmu[i,yr] <- x[(i-1),yr] + ifelse(CDDs[i,yr]<CDDtrigger,0,b0 + (-1 * b0 * x[(i-1),yr]))
  xmu[i,yr] <- x[(i-1),yr] + ifelse(Tair[i,yr]<baseTemp,b0 + (-1 * b0 * x[(i-1),yr]),0)
  #xmu[i,yr] <- x[(i-1),yr] + ifelse(D[i,yr]<Dtrigger,b0 + (-1 * b0 * x[(i-1),yr]),0)
  #xmu[i,yr] <- x[(i-1),yr] + b0 + (-1 * b0 * x[(i-1),yr])
  x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
  #x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
  x[i,yr] <- min(max(x1[i,yr],0),0.999)
  
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a,x1.b)
  #CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  #CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
  
  #Dtrigger ~ dunif(Dtrigger.lower,Dtrigger.upper)
  
  b0 ~ dunif(b0_lower,b0_upper)
  #b2 <- 0
  #b3 <- 0
  #b4 <- 0
  #b5 <- 0
  # b2 ~ dunif(b2_lower,b2_upper)
  # b3 ~ dunif(b3_lower,b3_upper)
  # b4 ~ dunif(b4_lower,b4_upper) 
  # b5 ~ dunif(b5_lower,b5_upper)
  }
  
  "
}else if(vars=="logistic"){ #Only b1
  outputFileName <- paste(siteName,"_CDD",baseTempOrig,"_logistic_noTrigger_calibration_varBurn.RData",sep="")
  
  
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 1
  
  
  variableNames <- c("p.PC","p.proc","x","b1","b2","CDDtrigger")
  
  generalModel = "
  model {
  ### Data Models for complete years
  for(yr in 1:(N)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N)){
  for(i in 2:n){
  
  ##Calculate CDD
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
  
  #xmu[i,yr] <- x[(i-1),yr] + ifelse(CDDs[i,yr]>0,(b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2),0)  
  xmu[i,yr] <- x[(i-1),yr] + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)
  x1[i,yr] ~ dnorm(xmu[i,yr],p.proc)
  x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
  
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a,x1.b)
  CDDs[1,yr] <- 0 ##Assumption built in based off of cut-off days
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  #CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
  
  b1 ~ dunif(b1_lower,b1_upper)
  b2 ~ dunif(b2_lower,b2_upper)
  }
  
  "
}
else if(vars=="noCov_oneYear"){ #full 
  #outputFileName <- paste(siteName,"_CDD",baseTemp,"_b3b4b5_noTrigger_TairShifted_calibration_varBurn.R",sep="")
  outputFileName <- paste(siteName,"_noCov_oneYear_noTrigger_calibration_varBurn2.RData",sep="")
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 1#0
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 1
  dataFinal$yr <- 1
  
  variableNames <- c("p.PC","p.proc","x","b0","b1","b2")#,"CDDtrigger")
  #variableNames <- c("p.PC","p.proc","x","b0","b1","b3","b4","b5")#,"CDDtrigger")
  
  generalModel = "
  model {
  ### Data Models for complete years
  for(i in 1:n){
  p[i,1] ~ dnorm(x[i,1],p.PC)
  }
  
  #### Process Model
  for(i in 2:n){
  xmu[i,yr] <- x[(i-1),yr] + b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)
  x[i,yr] ~ dnorm(xmu[i,yr],p.proc) T(0,1)
  #x[i,yr] <- min(max(x1[i,yr],0),x[(i-1),yr])
  }
  
  #### Priors
  x[1,yr] ~ dbeta(x1.a,x1.b)
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  
  b0 ~ dunif(b0_lower,b0_upper)
  b1 ~ dunif(b1_lower,b1_upper)
  b2 ~ dunif(b2_lower,b2_upper)
  }
  "
  
}


# library(PhenoForecast)
# library(PhenologyBayesModeling)
# library(rjags)
# library(runjags)
# library(suncalc)
# library(rnoaa)
# library(doParallel)
# 
# siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
# endDate <- as.Date("2019-12-31")
# dataDirectory <- "PhenologyForecastData/"
# 
# modelDescriptions <- data.frame(matrix(ncol=3,nrow=15))
# colnames(modelDescriptions) <- c("rowNum","baseTemp","vars")
# modelDescriptions$rowNum <- c(rep(1,3),rep(10,3),rep(13,3),rep(9,3),rep(6,3))
# modelDescriptions$baseTemp <- rep(c(10,15,20),5)
# #modelDescriptions$vars <- rep(c(rep("b3b4b5",3),rep("b3b4",3),rep("b3b5",3),rep("b4b5",3)),5)
# modelDescriptions$vars <- rep("hetero",15)
# 
# r=1
# i <- modelDescriptions$rowNum[r]
# baseTemp <- modelDescriptions$baseTemp[r]
# vars <- modelDescriptions$vars[r]
# siteName <- as.character(siteData[i,1])
# print(siteName)
# 
# lat <- as.numeric(siteData[i,2])
# long <- as.numeric(siteData[i,3])
# startDate <- (as.Date(siteData[i,7]))
# URL <- as.character(siteData$URL[i])
# URL2 <- as.character(siteData$URL2[i])
# URL3 <- as.character(siteData$URL3[i])
# if(nchar(URL2)>0){
#   URL <- c(URL,URL2)
#   if(nchar(URL3)>0){
#     URL <- c(URL,URL3)
#   }
# }
# TZ <- as.numeric(siteData[i,6])
# 
# rescaleFile <- paste(dataDirectory,siteName,"_forecast_phenoFits_PC.csv",sep="")
# rescaleData <- read.csv(rescaleFile,header=TRUE)
# cMeans.p <- rescaleData$cMeans.p
# dMeans.p <- rescaleData$dMeans.p
# #ERA5dataFolder <- paste("/projectnb/dietzelab/kiwheel/ERA5/Data/",siteName,"/",sep="")
# ERA5dataFolder <- dataDirectory
# URLs=URL
# cValsPC=cMeans.p
# dValsPC=dMeans.p
# 
# 
# ###Download PhenoCam data and format
# phenoData <- matrix(nrow=0,ncol=32)
# for(u in 1:length(URLs)){
#   print(URLs[u])
#   phenoDataSub <- download.phenocam(URLs[u])
#   phenoData <- rbind(phenoData,phenoDataSub)
# }
# ##Order and remove duplicate PC data
# phenoData2 <- phenoData[order(phenoData$date),]
# phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
# phenoData <- phenoData3
# 
# phenoData <- phenoData[phenoData$date<endDate,]
# p.old <- phenoData$gcc_mean
# time.old <-  as.Date(phenoData$date)
# 
# days <- seq(as.Date(startDate),(as.Date(endDate)),"day")
# p <- rep(NA,length(days))
# 
# for(i in 1:length(p.old)){
#   p[which(days==time.old[i])] <- p.old[i]
# }
# 
# months <- lubridate::month(days)
# years <- lubridate::year(days)
# 
# 
# #pdf(paste(siteName,"_",endDate,"_DataPlots.pdf",sep=""),width=10,height=8)
# 
# dat2 <- data.frame(dates=days,years=years,months=months,p=p)
# datTairEns <- load_ERA5_Tair_New(ERA5dataFolder=ERA5dataFolder)
# TairMu <- apply(X=datTairEns,MARGIN=2,FUN=mean)
# TairPrec <- 1/apply(X=datTairEns,MARGIN=2,FUN=var)
# dat2$TairMu <- TairMu
# dat2$TairPrec<- TairPrec
# 
# if(TZ==5){
#   TZ_name <- "America/New_York"
# }else if (TZ==6){
#   TZ_name <- "America/Chicago"
# }else{
#   print("Unrecognized TZ")
# }
# 
# dayLengths <- numeric()
# 
# for(d in 1:length(days)){
#   suntimes <- getSunlightTimes(date=days[d],
#                                lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),
#                                tz = TZ_name)
#   dayLengths <- c(dayLengths,as.numeric(suntimes$nauticalDusk-suntimes$nauticalDawn))
# }
# print("Day length processed")
# 
# dat2$D <- dayLengths
# 
# dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(213,365),]
# #dat2 <- dat2[lubridate::month(dat2$dates)%in%seq(7,12),]
# nrowNum <- 365-212
# p <- matrix(nrow=nrowNum,ncol=0)
# TairMu <- matrix(nrow=nrowNum,ncol=0)
# D <- matrix(nrow=nrowNum,ncol=0)
# 
# TairPrec <- matrix(nrow=nrowNum,ncol=0)
# valNum <- 0
# days2 <- matrix(nrow=nrowNum,ncol=0)
# for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
#   subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
#   valNum <- valNum + 1
# 
#   c <- cValsPC[valNum]
#   d <- dValsPC[valNum]
#   newCol <- PhenologyBayesModeling::rescale(yseq=subDat$p,c=c,d=d)
# 
#   p <- cbind(p,newCol)
#   days2 <- cbind(days2,as.Date(subDat$dates))
#   TairMu <- cbind(TairMu,subDat$TairMu)
#   #D <- cbind(D,subDat$D)
#   TairPrec <- cbind(TairPrec,subDat$TairPrec)
# 
# }
# p[p<0] <- 0
# p[p>1] <- 1
# #############
# 
# n <- nrowNum
# N <- ncol(p)
# CDDtrigger <- 0
# baseTemp <- 20
# Tair <- TairMu
# 
# b0 <- 0 #Probably unif(-1,0)
# b1 <- -1.2 #Probably unif(-1,1)
# b2 <- 0.2 #Probably unif(-1,1)
# b3 <- 0#Probably unif(0,1)
# b4 <- 0 #Probably unif(0,1)
# b5 <- 0 # Probably unif(0,1)
# 
# print("Data Processed")
# 
# yr <- 2
# CDDs <- matrix(ncol=N,nrow=n)
# x <- matrix(ncol=N,nrow=n)
# trigDeltaX <- matrix(ncol=N,nrow=n)
# deltaX <- matrix(ncol=N,nrow=n)
# 
# r2 <- 0.06
# 
# for(yr in 1:N){
#   CDDs[1,yr] <- 0
#   x[1,yr] <- 0.99
#   for(i in 2:n){
#     if(Tair[i,yr]<baseTemp){
#       CDDs[i,yr] <- CDDs[(i-1),yr] + baseTemp - Tair[i,yr]
#     }else{
#       CDDs[i,yr] <- CDDs[(i-1),yr]
#     }
#     #trigDeltaX[i,yr] <- max(0,x[(i-1),yr] - r2 * (1 - x[(i-1),yr]))
#     trigDeltaX[i,yr] <- b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2) + (b3 * Tair[i,yr]) + (b4 * D[i,yr]) + (b5 * Tair[i,yr] * D[i,yr])
#     if(CDDs[i,yr]>CDDtrigger){
#       deltaX[i,yr] <- trigDeltaX[i,yr]
#     }else{
#       deltaX[i,yr] <- 0
#     }
# 
#     x[i,yr] <- x[(i-1),yr] + deltaX[i,yr]
#     x[i,yr] <- min(max(0,x[(i-1),yr] + deltaX[i,yr]),1)
# 
#   }
# }
# #plot(seq(1,length(x[,yr])),x[,yr],typ="l")
# plot(seq(1,length(p)),p,pch=20)
# lines(seq(1,length(x)),x,col="green")
# 
# par(mfrow=c(2,1),mai=c(0.1,0.1,0.1,0.1))
# 
# par(mfrow=c(1,1),mai=c(0.1,0.1,0.1,0.1))
# 
# plot(seq(1,length(CDDs)),CDDs,pch=20)
# abline(h=CDDtrigger,col="red")



#####################
sq <- c(0.2,0.3,0.38,0.45,0.46,0.45,0.38,0.3,0.2)

sq1 <- c(0.4,0.3,0.22,0.18,0.178,0.18,0.22,0.3,0.4)
sq2 <- sq1*1.5
sq2[c(1,9)] <- 0.7
sq2[c(2,8)] <- 0.5
sq2 <- c(0.4,0.3,0.22,0.18,0.178,0.18,0.22,0.3,0.4)
plot(1:length(sq1),sq2,pch=20,ylim=c(0,1))
points(1:length(sq1),sq1,pch=20,col="cyan")
plot(1:length(sq1),(sq2-sq1)/(sq1+sq2),pch=20)

#sq2 <- sq*0.37-0.1
sq2 <- (sq+0.1)/0.37

sq3 <- sq-0.1
plot(1:length(sq1),sq1,pch=20)




siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)

days <- seq(1,365)
for(i in 1:length(days)){
  if(as.numeric(days[i])<10){
    days[i] <- paste("00",as.character(days[i]),sep="")
  }
  else if(as.numeric(days[i])<100){
    days[i] <- paste("0",days[i],sep="")
  }
}
year <- 2018

for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[s])
  
  emptyFiles <- 0
  for(day in days){
    #fileName <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/",siteName,"_GOES_diurnal_",year,day,".csv",sep="")
    fileName <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="")
    if(file.exists(fileName)){
      dayData <- read.csv(fileName,header=FALSE)
      if(ncol(dayData)<2){
        emptyFiles <- emptyFiles + 1
      }
    }
  }
  for(day in days){
    fileName <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/",siteName,"_GOES_diurnal_",year,day,".csv",sep="")
    #fileName <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/GOES_NDVI_Diurnal",siteName,"_",year,day,".csv",sep="")
    if(file.exists((fileName))){
      dayData <- read.csv(fileName,header=FALSE)
      if(ncol(dayData)<2){
        emptyFiles <- emptyFiles + 1
      }
    }
  }
  print(siteName)
  print(emptyFiles)
}

#############################

#!/usr/bin/env Rscript

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
library("ncdf4")
library(plyr)
library("PhenologyBayesModeling")
library(doParallel)
library("rjags")
library("runjags")


siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
savePath <- paste(getwd(),"/PhenologyForecastData/GOES_DiurnalFits/",sep="")
startDate=as.Date("2018-01-01")
endDate=as.Date("2018-12-31")
yr <- 2018

for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData[s,1])
  print(siteName)
  outDataFile <- paste(siteName,"_",yr,"_diurnalFitDataModel.RData",sep="")
  diurnalFits <- intersect(dir(path="PhenologyForecastData/GOES_DiurnalFits",pattern="varBurnFilter"),
                           dir(path="PhenologyForecastData/GOES_DiurnalFits",
                               pattern=paste(siteName,"_",yr,sep="")))
  c.vals <- numeric()
  prec.vals <- numeric()
  days <- numeric()
  cts <- numeric()
  for(i in 1:length(diurnalFits)){
    print(diurnalFits[i])
    load(paste("PhenologyForecastData/GOES_DiurnalFits/",diurnalFits[i],sep=""))
    if(typeof(var.burn)!=typeof(FALSE)){
      out.mat <- data.frame(as.matrix(var.burn))
      print(colnames(out.mat))
      c <- mean(out.mat$c)
      prec <- 1/var(out.mat$c)
      
      dy <- strsplit(diurnalFits[i],"_")[[1]][2]
      
      dayDataFile <- intersect(dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=siteName),
                               dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=paste(dy,".csv",sep="")))
      #dayDataFile <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/","GOES_NDVI_Diurnal",siteName,dy,".csv",sep="")
      print(dayDataFile)
      dayData <- read.csv(paste("PhenologyForecastData/GOES_NDVI_DiurnalData/",dayDataFile,sep=""),header=FALSE)
      ct <- length(dayData[2,][!is.na(dayData[2,])])
      if(ct>1){
        c.vals <- c(c.vals,c)
        prec.vals <- c(prec.vals,prec)
        days <- c(days,substr(dy,5,7))
        cts <- c(cts,ct)
      }
    }
    data <- list()
    data$x <- as.numeric(days)
    data$y <- as.numeric(c.vals)
    data$obs.prec <- as.numeric(prec.vals)
    data$n <- length(data$x)
    save(data,file=outDataFile)
    save(file=paste(siteName,"_",yr,"_counts.RData",sep=""),cts)
    print("Done with creating Data")
  }
  
}



#############################



#!/usr/bin/env Rscript

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
library("ncdf4")
library(plyr)
library("PhenologyBayesModeling")
library(doParallel)
library("rjags")
library("runjags")

#siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
savePath <- paste(getwd(),"/PhenologyForecastData/GOES_DiurnalFits/",sep="")
startDate=as.Date("2018-01-01")
endDate=as.Date("2018-12-31")
yr <- 2018

for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData[s,1])
  print(siteName)
  diurnalFits <- intersect(dir(path="PhenologyForecastData/GOES_DiurnalFits",pattern="varBurn"),
                           dir(path="PhenologyForecastData/GOES_DiurnalFits",
                               pattern=paste(siteName,"_",yr,sep="")))
  cts <- numeric()
  for(i in 1:length(diurnalFits)){
    #print(diurnalFits[i])
    load(paste("PhenologyForecastData/GOES_DiurnalFits/",diurnalFits[i],sep=""))
    if(typeof(var.burn)!=typeof(FALSE)){
      dy <- strsplit(diurnalFits[i],"_")[[1]][2]
      
      dayDataFile <- intersect(dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=siteName),
                               dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=paste(dy,".csv",sep="")))
      #dayDataFile <- paste("PhenologyForecastData/GOES_NDVI_DiurnalData/","GOES_NDVI_Diurnal",siteName,dy,".csv",sep="")
      #print(dayDataFile)
      dayData <- read.csv(paste("PhenologyForecastData/GOES_NDVI_DiurnalData/",dayDataFile,sep=""),header=FALSE)
      ct <- length(dayData[2,][!is.na(dayData[2,])])
      #print(ct)
      cts <- c(cts,ct)
    }
  }
}




##############################



# #!/usr/bin/env Rscript
# 
# #install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
# library("ncdf4")
# library(plyr)
# library("PhenologyBayesModeling")
# library(doParallel)
# library("rjags")
# library("runjags")
# 
# #siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
# siteData <- read.csv("PhenologyForecastData/GOES_Paper_Sites_FINAL.csv",header=TRUE)
# savePath <- paste(getwd(),"/PhenologyForecastData/GOES_DiurnalFits/",sep="")
# startDate=as.Date("2018-01-01")
# endDate=as.Date("2018-12-31")
# yr <- 2018
# 
# for(s in 2:nrow(siteData)){
#   siteName <- as.character(siteData[s,1])
#   print(siteName)
#   diurnalFits <- intersect(dir(path="PhenologyForecastData/GOES_DiurnalFits",pattern="varBurn"),
#                            dir(path="PhenologyForecastData/GOES_DiurnalFits",
#                                pattern=paste(siteName,"_",yr,sep="")))
#   cts <- numeric()
#   for(i in 1:length(diurnalFits)){
#     #print(diurnalFits[i])
#     load(paste("PhenologyForecastData/GOES_DiurnalFits/",diurnalFits[i],sep=""))
#     if(typeof(var.burn)!=typeof(FALSE)){
#       dy <- strsplit(diurnalFits[i],"_")[[1]][2]
#       
#       dayDataFile <- intersect(dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=siteName),
#                                dir(path="PhenologyForecastData/GOES_NDVI_DiurnalData",pattern=paste(dy,".csv",sep="")))
# 
#       dayData <- read.csv(paste("PhenologyForecastData/GOES_NDVI_DiurnalData/",dayDataFile,sep=""),header=FALSE)
#       ct <- length(dayData[2,][!is.na(dayData[2,])])
#       #print(ct)
#       if(ct>1){
#         cts <- c(cts,ct)
#       }
#     }
#   }
#   save(file=paste(siteName,"_",yr,"_counts.RData",sep=""),cts)
# }









library(devtools)

install_github("k-wheeler/NEFI_pheno/GOESDiurnalNDVI",build_vignettes = TRUE,force=TRUE)
library(GOESDiurnalNDVI)
browseVignettes("GOESDiurnalNDVI")

load("GOESDiurnalNDVI/RussellSage_2017234_varBurn.RData")
var.burn <- varBurn
save(var.burn,file="GOESDiurnalNDVI/RussellSage_2017234_varBurn.RData")

help(package = "GOESDiurnalNDVI", help_type = "html")
browseVignettes("GOESDiurnalNDVI")
install.packages("GOESDiurnalNDVI_0.1.0.tar.gz",repos=NULL,type="source")
library(GOESDiurnalNDVI)


##################################
library(rjags)
library(runjags)
library("PhenoForecast")
library("PhenologyBayesModeling")
library(ecoforecastR)

siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
forecastDataFolder <- "PhenologyForecastData/ForecastOutputs/AllForecasts/"
dataDirectory <- "PhenologyForecastData/"
Nmc <- 10000 #Number of model runs
#allDates <- seq(as.Date("2019-01-23"),as.Date("2019-06-06"),"day")
allDates <- c(seq(as.Date("2019-01-23"),as.Date("2019-01-25"),"day"),
              as.Date("2019-02-03"),as.Date("2019-02-05"),
              seq(as.Date("2019-02-07"),as.Date("2019-06-06"),"day"))
#allDates <- allDates[1:(length(allDates)/2)]
#allDates <- allDates[(length(allDates)/2):length(allDates)]
i <- 6

##General site-specific info
siteName <- as.character(siteData$siteName[i])
print(siteName)
plotFolder <- paste("VPplots/",siteName,"/",sep="")
URL <- as.character(siteData$URL[i])
lat <- as.numeric(siteData$Lat[i])
long <- as.numeric(siteData$Long[i])
station <- as.character(siteData$metStation[i])
siteStartDate <- as.character(siteData$startDate[i])
#pdf(file=paste(siteName,"_PhenologyForecast_variancePartition.pdf",sep=""),height=10,width=6)
#par(mfrow=c(3,1))

allDates <- c(as.Date("2019-04-07"),
              as.Date("2019-04-14"),
              as.Date("2019-04-21"),
              as.Date("2019-04-23"),
              as.Date("2019-04-26"))
##Date info 
for(d in 1:length(allDates)){
  
  calEndDate <- allDates[d]
  #calEndDate <- as.Date("2019-01-23")#as.Date("2019-05-20")
  #calEndDate <- as.Date("2019-05-20")
  print(calEndDate)
  forStartDate <- calEndDate + 1
  forEndDate <- (forStartDate+13)
  forDates <- seq(forStartDate,forEndDate,"day")
  yrDates <- seq(as.Date("2019-01-01"),forEndDate,"day")
  
  forecastLength <- 14
  plotDates <- as.numeric(format(forDates,"%j"))
  
  LCfileName <- paste(forecastDataFolder,siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
  #if(file.exists(LCfileName)){
  ##Assemble Sf forecast
  #GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/",siteName,"/",calEndDate,"/",sep="")
  #GEFS_Files <- dir(path=GEFS_Directory,pattern="NOAA_GEFS")
  
  #SfsALL <- createSf(lat=lat,long=long,dates=yrDates,siteName=siteName,
  #                   dataDirectory=dataDirectory,endDate=forEndDate,
  #                   GEFS_Files=GEFS_Files,GEFS_Directory=GEFS_Directory,
  #                   forecastLength=forecastLength, station=station,
  #                   calDatesT = FALSE) 
  
  #save(SfsALL,file=paste("SfsALL_",calEndDate,".RData",sep=""))
  load(paste("SfsALL_",calEndDate,".RData",sep=""))
  #load("SfsALL_140.RData")
  #load("SfsALL_140_new.RData")
  #load("harvard_2009-01-01_2019-05-20_LC2_outBurn.RData")
  #outBurnLC2 <- outBurnLC
  #Sf.means <- matrix(colMeans(SfsALL),1,forecastLength)
  
  #LCfileName <- paste(siteName,"_",siteStartDate,"_",calEndDate,"_LC2_outBurn.RData",sep="")
  
  load(LCfileName)
  out.mat.par <- data.frame(as.matrix(outBurnLC2$params))
  dayNumber <- dim(as.matrix(outBurnLC2$predict))[2]-14
  out.mat.LC <- as.matrix(outBurnLC2$predict)
  out.mat.LC <- out.mat.LC[,1:dayNumber]
  ICs <- out.mat.LC[,dayNumber]
  
  ##Deterministic prediction
  LC.det <- forecastLogCov2(IC=mean(ICs),
                            trans = mean(out.mat.par$trans),
                            b1 = mean(out.mat.par$b1),
                            Sf = matrix(SfsALL$Sf[(length(SfsALL$Sf)-13):length(SfsALL$Sf)],1,14),
                            Q=0,
                            n=1,
                            NT=14)
  
  ##Initial Condition Uncertainty
  rndNums <- sample.int(nrow(out.mat.par),Nmc,replace=TRUE)
  LC.I <- forecastLogCov2(IC=ICs[rndNums],
                          trans = mean(out.mat.par$trans),
                          b1 = mean(out.mat.par$b1),
                          Sf = matrix(SfsALL$Sf[(length(SfsALL$Sf)-13):length(SfsALL$Sf)],1,14),
                          Q=0,
                          n=Nmc,
                          NT=14)
  LC.I.ci <- apply(LC.I,2,quantile,c(0.025,0.5,0.975))
  
  ##Parameter Uncertainty 
  LC.IP <- forecastLogCov2(IC=ICs[rndNums],
                           trans = out.mat.par$trans[rndNums],
                           b1 = out.mat.par$b1[rndNums],
                           Sf = matrix(SfsALL$Sf[(length(SfsALL$Sf)-13):length(SfsALL$Sf)],1,14),
                           Q=0,
                           n=Nmc,
                           NT=14)
  LC.IP.ci = apply(LC.IP,2,quantile,c(0.025,0.5,0.975))
  
  ##Driver Uncertainty
  SfsSamp <- matrix(nrow=Nmc,ncol=0)
  
  for(s in (length(SfsALL$Sf)-13):length(SfsALL$Sf)){ ##Only need the last 14
    SfsSamp <- cbind(SfsSamp,rnorm(Nmc,SfsALL$Sf[s],1/sqrt(SfsALL$Sfprec[s])))
  }
  LC.IPD <- forecastLogCov2(IC=ICs[rndNums],
                            trans = out.mat.par$trans[rndNums],
                            b1 = out.mat.par$b1[rndNums],
                            Sf = SfsSamp,
                            Q=0,
                            n=Nmc,
                            NT=14)
  LC.IPD.ci = apply(LC.IPD,2,quantile,c(0.025,0.5,0.975))
  
  ##Process Uncertainty
  Qmc <- 1/sqrt(out.mat.par[rndNums,"p.proc"])
  LC.IPDE <- forecastLogCov2(IC=ICs[rndNums],
                             trans = out.mat.par$trans[rndNums],
                             b1 = out.mat.par$b1[rndNums],
                             Sf = SfsSamp,
                             Q=Qmc,
                             n=Nmc,
                             NT=14)
  LC.IPDE.ci <- apply(LC.IPDE,2,quantile,c(0.025,0.5,0.975))
  
  varI     <- apply(LC.I,2,var)
  varIP    <- apply(LC.IP,2,var)
  varIPD   <- apply(LC.IPD,2,var)
  varIPDE  <- apply(LC.IPDE,2,var)
  varMat   <- rbind(varI,varIP,varIPD,varIPDE)
  
  ## out-of-sample stacked area plot
  V.pred.rel <- apply(varMat,2,function(x) {x/max(x)})
  
  #V.pred.rel <- matrix(ncol=ncol(varMat),nrow=0)
  #for(c in 1:ncol(varMat)){
  #  x <- varMat[,c]
  #  V.pred.rel <- rbind(V.pred.rel,x/max(x))
  
  #}
  #plotFileName <- paste(plotFolder,"LC2/",siteName,"_LC2_",calEndDate,"_varPartition.png",sep="")
  #png(file=plotFileName, width=10, height=5,units="in",res=1000)
  plot(forDates,V.pred.rel[1,],ylim=c(0,1),type='n',main="Relative Variance Logistic with Covariate",ylab="Proportion of Variance",xlab="time")
  ciEnvelope(forDates,V.pred.rel[3,],V.pred.rel[4,],col="blue")
  ciEnvelope(forDates,rep(0,ncol(V.pred.rel)),V.pred.rel[1,],col="gray")
  ciEnvelope(forDates,V.pred.rel[1,],V.pred.rel[2,],col="green")
  ciEnvelope(forDates,V.pred.rel[2,],V.pred.rel[3,],col="pink")
  
  #dev.off()
}








########Investigating latent phenologicals states in spring between different sources
startDay <- 182
endDay <- 365+181

xseq <- seq(366,endDay)
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
#pdf("TransFigure.pdf",width=60,height=30)
#par(mfrow=c(2,3))
i <- 1
iseq <- c(seq(1,6),seq(8,10),seq(15,20))

library("PhenologyBayesModeling")

plot(numeric(),numeric(),ylim=c(0,1),xlim=c(0,1),xlab="PhenoCam",ylab="GOES")
for(i in iseq){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  
  ##GOES
  inputFileName <- paste(siteName,"_overall6_varBurn.RData",sep="")
  load(inputFileName)
  
  #var.mat<-as.matrix(GOES.md.out)
  var.mat <- data.frame(as.matrix(var.burn))
  rndNums <- sample(1:nrow(var.mat),10000,replace=T)
  TranS <- var.mat$TranS[rndNums]
  bS <- var.mat$bS[rndNums]
  c <- var.mat$c[rndNums]
  d <- var.mat$d[rndNums]
  k <- var.mat$k[rndNums]
  ycred.GOES <- matrix(0,nrow=10000,ncol=length(xseq))
  for(g in 1:10000){
    Ey <- pheno.logistic(Tran=TranS[g],b=bS[g],c=c[g],d=d[g],xseq=xseq)
    Ey <- rescale(c=c[g],d=d[g],yseq = Ey)
    ycred.GOES[g,] <- Ey
  }
  
  ##MODIS NDVI
  inputFileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_NDVI_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-data.frame(as.matrix(MODIS.N.md.out))
  
  rndNums <- sample(1:nrow(var.mat),10000,replace=T)
  TranS <- var.mat$TranS[rndNums]
  bS <- var.mat$bS[rndNums]
  c <- var.mat$c[rndNums]
  d <- var.mat$d[rndNums]
  k <- var.mat$k[rndNums]
  ycred.MODISN <- matrix(0,nrow=10000,ncol=length(xseq))
  for(g in 1:10000){
    Ey <- pheno.logistic(Tran=TranS[g],b=bS[g],c=c[g],d=d[g],xseq=xseq)
    Ey <- rescale(c=c[g],d=d[g],yseq = Ey)
    ycred.MODISN[g,] <- Ey
  }
  
  ##MODIS EVI
  inputFileName <- paste(siteName,"_",startDay,"_",endDay,"_MODIS_DQF_EVI_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-data.frame(as.matrix(MODIS.E.md.out))
  
  rndNums <- sample(1:nrow(var.mat),10000,replace=T)
  TranS <- var.mat$TranS[rndNums]
  bS <- var.mat$bS[rndNums]
  c <- var.mat$c[rndNums]
  d <- var.mat$d[rndNums]
  k <- var.mat$k[rndNums]
  ycred.MODISE <- matrix(0,nrow=10000,ncol=length(xseq))
  for(g in 1:10000){
    Ey <- pheno.logistic(Tran=TranS[g],b=bS[g],c=c[g],d=d[g],xseq=xseq)
    Ey <- rescale(c=c[g],d=d[g],yseq = Ey)
    ycred.MODISE[g,] <- Ey
  }
  
  ##PhenoCam
  inputFileName <- paste(siteName,"_",startDay,"_",endDay,"_PC_varBurn.RData",sep="")
  load(inputFileName)
  
  var.mat<-data.frame(as.matrix(PC.md.out))
  rndNums <- sample(1:nrow(var.mat),10000,replace=T)
  TranS <- var.mat$TranS[rndNums]
  bS <- var.mat$bS[rndNums]
  c <- var.mat$c[rndNums]
  d <- var.mat$d[rndNums]
  k <- var.mat$k[rndNums]
  ycred.PC <- matrix(0,nrow=10000,ncol=length(xseq))
  for(g in 1:10000){
    Ey <- pheno.logistic(Tran=TranS[g],b=bS[g],c=c[g],d=d[g],xseq=xseq)
    Ey <- rescale(c=c[g],d=d[g],yseq = Ey)
    ycred.PC[g,] <- Ey
  }
  points(apply(ycred.PC,2,quantile,0.50),apply(ycred.GOES,2,quantile,0.50),pch=20)
  #plot(xseq,apply(ycred.PC,2,quantile,0.50),pch=20,col=PC.col)
  #points(xseq,apply(ycred.GOES,2,quantile,0.50),pch=20,col=GOES.col)
  #abline(0,1,col="red")
  
}




#################################
library("PhenoForecast")

library("dplyr")


siteData <- read.csv("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
dataDirectory="/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyForecastData/"


##Need to make plots of the different met sources versus each other:
##ERA5 versus NOAA met

##NOAA met vs GEFS Forecast

##ERA5 versus willowCreek met

##willowCreek met vs NOAA met

##GEFS versus willowCreek

iseq <- c(seq(1,6),8,9,10,15,16)
endDate <- as.Date("2019-04-04")
pdf("PhenologyForecast_Met_ComparisonPlots.pdf",height=6,width=6)
par(pty="s")
for(i in iseq){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  lat <- as.numeric(siteData$Lat[i])
  long <- as.numeric(siteData$Long[i])
  startDate <- as.Date(siteData$startDate[i])
  station <- as.character(siteData$metStation[i])
  calDates <- seq(startDate,as.Date("2018-12-31"),"day")
  
  GEFS_Directory <- paste("/projectnb/dietzelab/WeatherForecast/NOAA_GEFS/Data/",siteName,"/",endDate,"/",sep="")
  GEFS_Files <- dir(path=GEFS_Directory,pattern="NOAA_GEFS")
  
  TairsERA5 <- load_ERA5_Tair(lat=lat,long=long,years=seq(lubridate::year(startDate),2018)) ##columns are each an ensemble (not divided by year)
  TairsNOAAMet <- load_NOAA_met(station=station,startDate=startDate,endDate=endDate) ##Array of numeric values
  
  TairsGEFS <- numeric()
  for(e in 1:length(GEFS_Files)){
    TairsForecastInd <- load_GEFS_Forecast(dataDirectory=GEFS_Directory,fileName=GEFS_Files[e])
    TairsGEFS <- cbind(TairsForecast,TairsForecastInd)
  }
  
  ##Need mean and sd for ERA5 and GEFS
  TairsERA5.mu <- rowMeans(TairsERA5)
  TairsERA5.sd <- sqrt(apply(TairsERA5,1,var))
  
  TairsGEFS.mu <- colMeans(TairsGEFS)
  TairsGEFS.sd <- sqrt(apply(TairsGEFS,2,var))
  mdl <- lm(TairsNOAAMet[1:length(TairsERA5.mu)]~TairsERA5.mu)
  
  plot(TairsERA5.mu,TairsNOAAMet[1:length(TairsERA5.mu)],pch=20,main=paste(siteName,"NOAA met versus ERA5"))
  abline(1,1,col="red")
  abline(mdl,col="cyan")
  
  # if(siteName=="willowcreek"){
  #    TairsCalIndWC <- download_US_WCr_met(start_date=calDates[1],end_date=calDates[length(calDates)])
  #    TairsCurrentIndWC <- download_US_WCr_met(start_date=as.Date("2019-01-01"),end_date=(endDate-forecastLength))
  # }
  plot(seq(1,length(TairsERA5.mu)),TairsERA5.mu,pch=20,main=siteName)
  points(seq(1,length(TairsNOAAMet[1:length(TairsERA5.mu)])),TairsNOAAMet[1:length(TairsERA5.mu)],pch=20,col="red")
  
  
}
dev.off()
############################







#
NDSI.files <- intersect(dir(pattern="NDSI"),dir(pattern="shenandoah"))
output2 <- matrix(ncol=3,nrow=0)
for(i in 1:length(NDSI.files)){
  dat <- read.csv(NDSI.files[i],header=FALSE)
  #print(length(na.omit(as.numeric(dat[2,]))))
  if(length(na.omit(as.numeric(dat[2,])))>0){
    day <- strsplit(NDSI.files[i],split="_")[[1]][4]
    #print(NDSI.files[i])
    subMean <- mean(na.omit(as.numeric(dat[2,])))
    output2 <- rbind(output2,c(day,subMean,length(na.omit(as.numeric(dat[2,])))))
  }
}
output2

##Investigating hdfeos files for the NWP snow masks
install.packages("ncdf4",configure.args = --enable-hdf4)
library("ncdf4")
library("ncdf4.helpers")

nc<- nc_open("NISE_SSMISF18_20180101.HDFEOS")

v1 <- nc$var[["Extent"]]
z_all <- ncvar_get(nc, v1)
zv <- as.vector(as.single(z_all))
zz <- file("tmpbin", "wb")
writeBin(zv, zz)
close(zz)
zz2 <- file("tmpbin", "rb")
#zs <- readBin(zz2, numeric(), length(zv), endian="little")

zs <- readBin(zz2, numeric(), size=4, length(zv), endian="little")
zs <- readBin(con=zz2, what=numeric(), size=4, n=length(zv), endian="little")
zs <- readBin(con=zz2, what=numeric(), n=length(zv), endian="little",signed=FALSE)

close(zz2)
dim(zs) <- dim(z_all)
range(zs,na.rm=TRUE)
zs1 <- zs[zs>-Inf]
zs1 <- zs1[zs1<Inf]


##Need to investigate the Sf values to see how they increase
##(Plot how it changes throughout the year)

##Think of how it might compare to the timings of when you might expect the rate of increase to start increasing
###This should be related to the threshold value 
##Think about and research the best functional relationship for how the two relate (Need to research this relationship more)

##Probably could plot the %canopy values vs the Sf to get a better understanding of potential values for the ICs

##Need to compute an expected prior distribution based off of F.star values given in Melaas et al.2016
sp1 <- rnorm(10000,mean=171.5,sd=(193.7-153.6)/3.92)
sp2 <- rnorm(10000,mean=135.7,sd=(148.1-126.3)/3.92)
sp3 <- rnorm(10000,mean=261.6,sd=(290.1-237.2)/3.92)
sp4 <- rnorm(10000,mean=150.1,sd=(158.1-139.4)/3.92)
sp5 <- rnorm(10000,mean=165.6,sd=(185.1-152)/3.92)
mean(c(sp1,sp2,sp3,sp4,sp5))
sqrt(var(c(sp1,sp2,sp3,sp4,sp5)))

x <- seq(1,200,1)
r=1.8
#color <- rep(NA,200)
x <- rep(NA,200)
x[1] <- 0.01
for(i in 2:200){
  #print(i)
  x[i] <- x[i-1]+r*x[i-1]*(1-x[i-1])
}
plot(seq(1,200,1),x,pch=20)
abline(v=7,col="cyan")

library("ncdf4")
dataDirectory <- "phenologyForecastData/"
fileName <- "NOAA_GEFS.willowcreek.1.2019-01-25T12:00.2019-02-10T12:00.nc"
startDateTime <- strsplit(fileName,split="[.]")[[1]][4]
startTime <- strsplit(startDateTime,"T")[[1]][2]
startDate <-strsplit(startDateTime,"T")[[1]][1] 
weatherFile <- nc_open(paste(dataDirectory,fileName,sep=""))
Tair <- ncvar_get(weatherFile,"air_temperature")
timeVals <- weatherFile$dim$time$vals
dates <- as.POSIXlt(paste(startDate,startTime,sep=" "))+(timeVals*3600)
dates
currentDate <- lubridate::date(dates[1])

Tairs <- numeric()
subTair <- numeric()
for(i in 1:length(dates)){
  if(lubridate::date(dates[i])==currentDate){
    subTair <- c(subTair,Tair[i])
  }
  else{
    Tairs <- c(Tairs,mean(subTair))
    currentDate <- lubridate::date(dates[i])
    subTair <- Tair[i]
  }
  
}
Tairs <- c(Tairs,mean(subTair))
if(startTime!="00:00"){ ##Removes the partial start and end dates
  Tairs <- Tairs[2:(length(Tairs)-1)]
}


gefsFile <- nc_open("PhenologyForecastData/NOAA_GEFS.willowcreek.1.2019-01-25T12:00.2019-02-10T12:00.nc")
gfs <- structure(gefsFile)
length(ncvar_get(gefsFile,"air_temperature"))
ncdim_get(gefsFile,"time")

###########################################
library("PhenologyBayesModeling")
library("MODISTools")
siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
for(i in 1:nrow(siteData)){
  siteName <- as.character(siteData[i,1])
  print(siteName)
  lat <- as.character(siteData[i,2])
  long <- as.character(siteData[i,3])
  PFT <- as.character(siteData[i,5])
  if(PFT=="DB"){
    startDay <- 182
    endDay <- 365+181
  }
  else if(PFT=="SH"){
    startDay <- 110
    endDay <-90+365
  }
  MODIS_data(siteName = siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="EVI")
}




###########
library("phenocamr")
download_phenocam(site="harvard$",veg_type = "DB",out_dir="PhenoCam_Data")
list_sites()

#################
library("ncdf4")
library("ncdf4.helpers")
library(plyr)
#install.packages("ncdf4")

#fileName <- "OR_ABI-L2-ACMC-M3_G16_s20180011727198_e20180011729571_c20180011730171.nc"
fileName <- "GOES_Data/GOES_Data2017/OR_ABI-L1b-RadC-M3C05_G16_s20173001502209_e20173001504582_c20173001505029.nc"
fileName <- "GOES_Data2017/OR_ABI-L1b-RadC-M3C05_G16_s20180011502198_e20180011504571_c20180011505024.nc"
fileName <- "GOES_Data/GOES_Data2017/OR_ABI-L2-ACMC-M3_G16_s20173651952196_e20173651954569_c20173651955185.nc"
fileName <- 'GOES_Data/GOES_Data2017/OR_ABI-L1b-RadC-M3C03_G16_s20171821922189_e20171821924564_c20171821925010.nc'
fileName <- 'GOES_Data/GOES_Data2017/OR_ABI-L1b-RadC-M3C02_G16_s20171821922189_e20171821924562_c20171821924596.nc'
#fileName <- "NISE_SSMISF18_20180105.hdfeos"
#fileName <- "NISE_SSMISF18_20180105.HDFEOS"

C03.file <-nc_open(fileName)


filestrC05 <- paste("OR_ABI-L1b-RadC-M3C05_G16_s",day.time,sep="")
R5 <- ncvar_get(C05.file,"Rad")
R5.DQF <- ncvar_get(C05.file,"DQF")

structure(C05.file)

test <- h5file(ACM.file, 'r')

clouds <- ncvar_get(ACM.file,"BCM")[ACM.ind[1],ACM.ind[2]] #BCM= binary cloud mask
clouds.DQF <- ncvar_get(ACM.file,"DQF")[ACM.ind[1],ACM.ind[2]]
snowMask <- ncvar_get(ACM.file,"input_dynamic_ancillary_NWP_snow_mask_data")

################
diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  #print(k)
  bk <- which(round(xseq,digits=1)==k)
  #print(bk)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  #print(length(c(left,right)))
  return(c(left,right))
}

diurnalFiles <- dir(path="dailyNDVI_GOES",pattern = "Leftover")
length(diurnalFiles)
for(i in 1:length(diurnalFiles)){
  #print(i)
  dat <- read.csv(paste("dailyNDVI_GOES/",diurnalFiles[i],sep=""),header=FALSE)
  #print(dim(dat))
  if(is.na(dat[1,1])){
    print(diurnalFiles[i])
  }
}

#!/usr/bin/env Rscript

install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("MODISTools",repo="https://cloud.r-project.org/")
#install.packages("doParallel",repo="https://cloud.r-project.org/")
library("PhenologyBayesModeling")
library("rjags")
library("runjags")

library("ecoforecastR")


diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  #print(k)
  bk <- which(round(xseq,digits=1)==k)
  #print(bk)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  #print(length(c(left,right)))
  return(c(left,right))
}

siteName <- "russellSage"
xseq <- seq(0,25,0.1)

#outputFileName <- paste(siteName,"_diurnalExamples2.pdf",sep="")
#pdf(file=outputFileName,width=10,height=5)
#par(mfrow=c(1,2),mai=c(1,0.4,1,0.1))
#day.seq <- c("186","189","201","212","231","245","251","278","290","012","013","030","074","091","167","168")
#day.seq <- c("189","201","231","245","251","012","013","030")
#day.seq <- c("186","201","245","251","278","030","012","013","290")
#dy="030"
dy <- 231
#fileName <- dir(path="dailyNDVI_GOES",pattern=paste(dy,"_varBurn.RData",sep=""))
fileName <- intersect(dir(pattern=paste(dy,"_varBurn2.RData",sep="")),dir(pattern=siteName))
print(fileName)
load(fileName)
if(as.numeric(dy)<182){
  yr <- "2018"
}
else{
  yr <- "2017"
}
dat <- read.csv(paste("dailyNDVI_GOES/GOES_Diurnal_",siteName,"_",yr,dy,".csv",sep=""),header=FALSE)
out.mat <- as.matrix(var.burn)
a <- out.mat[,1]
rndNums <- sample(1:length(a),10000,replace=T)
a <- a[rndNums]
c <- out.mat[rndNums,2]
k <- out.mat[rndNums,3]
ycred <- matrix(0,nrow=10000,ncol=length(xseq))
for(g in 1:10000){
  Ey <- diurnalExp(a=a[g],c=c[g],k=k[g],xseq=xseq)
  ycred[g,] <- Ey
}
ci <- apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE)
date <- as.Date(as.numeric(dy),origin=as.Date(paste(as.character(as.numeric(yr)-1),"-12-31",sep="")))
plot(x=list(),y=list(),ylim=c(0,1),xlim=c(5,22),ylab="NDVI",xlab="Time (Hours)",cex.axis=2,bty='n',cex.lab=2)#,xaxt='n',yaxt='n')
#polygon(x=c(10,14,14,10),y=c(-1,-1,1.2,1.2),col="lightgray",border=NA)
#ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
lines(xseq,ci[2,],col="red",lwd=3)
points(as.numeric(dat[3,]),as.numeric(dat[2,]),pch=20)
abline(v=mean(k),lty=2)
#abline(v=quantile(k,0.025),lty=2)
#abline(v=quantile(k,0.975),lty=2)
abline(h=mean(c),lty=2)
abline(h=quantile(c,0.025),lty=2)
abline(h=quantile(c,0.975),lty=2)
#xseq2 <- seq(0,(mean(k),0.1)
lines(xseq,diurnalExp(a=quantile(a,0.025),k=mean(k),c=mean(c),xseq=xseq),lty=2)
lines(xseq,diurnalExp(a=quantile(a,0.975),k=mean(k),c=mean(c),xseq=xseq),lty=2)
}
#[127:251]
#dev.off()

siteName <- "russellSage"
diurnal.files <- intersect(dir(path="dailyNDVI_GOES",pattern=siteName),dir(path="dailyNDVI_GOES",pattern="Diurnal"))
days <- numeric()
maxs <- numeric()
i=1
for(i in 1:length(diurnal.files)){
  print(diurnal.files[i])
  dat <- read.csv(paste("dailyNDVI_GOES/",diurnal.files[i],sep=""),header=FALSE)
  NDVIs <- na.omit(as.numeric(dat[2,]))
  print(length(NDVIs))
  if(length(NDVIs)==0){
    maxs <- c(maxs,NA)
  }
  else{
    print(max(NDVIs))
    maxs <- c(maxs,max(NDVIs))
  }
  days <- c(days,substr(diurnal.files[i],30,32))
}
for(i in 1:length(days)){
  if(as.numeric(days[i])<182){
    days[i] <- as.numeric(days[i]) + 365
  }
}
plot(days,maxs)
final.dat <- rbind(days,maxs)
outFileName <- paste("GOES_NDVI_",siteName,"_",startDate,"_",endDate,"_Max.csv",sep="")
write.table(final.dat,file=outFileName,row.names=FALSE,col.names = FALSE,sep=",")

#install.packages("/projectnb/dietzelab/kiwheel/NEFI_pheno/PhenologyBayesModeling",repo=NULL)
#install.packages("MODISTools",repo="https://cloud.r-project.org/")
#install.packages("doParallel",repo="https://cloud.r-project.org/")
library("PhenologyBayesModeling")
library("rjags")
library("runjags")

library("ecoforecastR")
library("grDevices")


diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  #print(k)
  bk <- which(round(xseq,digits=1)==k)
  #print(bk)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  #print(length(c(left,right)))
  return(c(left,right))
}

siteName <- "russellSage"
xseq <- seq(0,25,0.1)

#outputFileName <- paste(siteName,"_diurnalExamples2.pdf",sep="")
#pdf(file=outputFileName,width=10,height=5)


par(mfrow=c(2,4),mai=c(0.3,0.3,0.4,0.05))
#day.seq <- c("186","189","201","212","231","245","251","278","290","012","013","030","074","091","167","168")
day.seq <- c("189","201","234","245","251","012","013","030")
#day.seq <- c("186","201","245","251","278","030","012","013","290")
#dy="030"
dy <- "234"
siteName <- "russellSage" 
#fileName <- dir(path="dailyNDVI_GOES",pattern=paste(dy,"_varBurn.RData",sep=""))
fileName <- intersect(dir(pattern=paste(dy,"_varBurn2.RData",sep="")),dir(pattern=siteName))
print(fileName)
load(fileName)
if(as.numeric(dy)<182){
  yr <- "2018"
}else{
  yr <- "2017"
}
dat <- read.csv(paste("dailyNDVI_GOES/GOES_Diurnal_",siteName,"_",yr,dy,".csv",sep=""),header=FALSE)
out.mat <- as.matrix(var.burn)

c <- out.mat[,2]
max.val <- max(na.omit(as.numeric(dat[2,])))
max.val

midday <- intersect(dat[,as.numeric(dat[3,])>10],dat[,as.numeric(dat[3,])<14])


date <- as.Date(as.numeric(dy),origin=as.Date(paste(as.character(as.numeric(yr)-1),"-12-31",sep="")))

#polygon(x=c(10,14,14,10),y=c(-1,-1,1.2,1.2),col="lightgray",border=NA)
#ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
#lines(xseq,ci[2,],col="black")

c.qnts <- quantile(c,c(0.025,0.975))
avg <- mean(as.numeric(midday[2,]),na.rm = TRUE)
avg.sd <- sqrt(var(as.numeric(midday[2,]),na.rm=TRUE))
avg.bot <- avg-1.96*(avg.sd/sqrt(length(na.omit(as.numeric(midday[2,])))))
avg.top <- avg+1.96*(avg.sd/sqrt(length(na.omit(as.numeric(midday[2,])))))
#avg.qnts <- quantile(midday[2,],c(0.025,0.975),na.rm=TRUE)
noon1 <- dat[2,119]
noon2 <- dat[2,120]
plot(x=list(),y=list(),main=date,ylim=c(0,0.9),xlim=c(5,20),ylab="NDVI",xlab="Time (Hours)",cex=2.5)
polygon(x=c(-1,-1,25,25),y=c(avg.bot,avg.top,avg.top,avg.bot),col=adjustcolor("#abd9e9",alpha.f=0.4),border=NA)
polygon(x=c(-1,-1,25,25),y=c(c.qnts[1],c.qnts[2],c.qnts[2],c.qnts[2]),col=adjustcolor("#d7191c",alpha.f=0.3),border=NA)

points(as.numeric(dat[3,]),as.numeric(dat[2,]),pch=20)
abline(h=avg,col="#abd9e9",lwd=2)

abline(h=mean(c),col="#d7191c")

abline(h=max.val,col="#fdae61")
"#2c7bb6"
abline(h=noon1,col="purple",lwd=1)
abline(h=noon2,col="purple",lwd=1)
points(as.numeric(dat[3,]),as.numeric(dat[2,]),pch=20)
legend("bottomright",col=c("#fdae61","purple","#d7191c","#abd9e9"),legend=c("Maximum","Potential Noons","Diurnal Fit","Window Average"),pch=c("-","-","-","-"),lwd=c(2,2,2,2))

#abline(v=12,col="red")

dev.off()





##########################
load("bullshoals_001_varBurn2.RData")
out.mat <- as.matrix(var.burn)
a <- out.mat[,1]
rndNums <-  sample(1:length(a),10000,replace=T)
a <- a[rndNums]
c <- out.mat[rndNums,2]
k <- out.mat[rndNums,3]
ycred <- matrix(0,nrow=10000,ncol=length(xseq))
for(g in 1:10000){
  Ey <- diurnalExp(a=a[g],c=c[g],k=k[g],xseq=xseq)
  ycred[g,] <- Ey
}
ci <- apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE)
year <- 2018
i <- "001"
siteName <- "bullshoals"
fileName <- paste("dailyNDVI_GOES/","GOES_Diurnal_",siteName,"_",year,i,".csv",sep="")
print(fileName)
dat <- read.csv(fileName,header=FALSE)
data <- list()
print(dim(dat))
data$x <- as.numeric(dat[3,])
data$y <- as.numeric(dat[2,])
plot(data$x,data$y,pch=20,ylim=c(0,1))
ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
points(data$x,data$y,pch=20)
abline(v=mean(k),col="red")
abline(v=dawnTime)
abline(v=duskTime)

library("ecoforecastR")
library("rjags")
library("runjags")

diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  #print(k)
  bk <- which(round(xseq,digits=1)==k)
  #print(bk)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  #print(length(c(left,right)))
  return(c(left,right))
}










#########################
library("PhenologyBayesModeling")
dat <- read.csv("GOES_NDVI_burns_2017-04-20_2018-03-31_Avg.csv",header=FALSE)
siteName <- "jerbajada"
#xseq <- seq(dat[1,1],dat[1,ncol(dat)],1)
#plot(as.numeric(dat[1,]),as.numeric(dat[2,]),pch=20)
#lines(xseq,shrublandYvals(Tran=114,b=-0.11,c=0.12,d=0.18,r=-0.018,k=150,xseq=xseq),col="red")

middayData <- read.csv("GOES_NDVI_jerbajada_2017-04-20_2018-03-31_Avg.csv",header=FALSE)
xseq <- seq(middayData[1,1],middayData[1,ncol(middayData)],1)
for(i in 1:length(middayData)){
  if(middayData[1,i]<110){
    middayData[1,i] <- middayData[1,i]+365
  }
}
middayData <- middayData[,!is.na(middayData[2,])]
middayData <- middayData[,!is.na(middayData[3,])]

middayDat <- matrix(ncol=ncol(middayData),nrow=5)
rownames(middayDat) <- c("Day","mean","sd","Q1","Q3")
middayDat[1,] <- as.numeric(middayData[1,])
middayDat[2,] <- as.numeric(middayData[2,])
middayDat[3,] <- sqrt(1/as.numeric(middayData[3,]))/(sqrt(as.numeric(middayData[4,])))
middayDat[4,] <- middayDat[2,]-middayDat[3,]
middayDat[5,] <- middayDat[2,]+middayDat[3,]

midday.dates <- as.Date(as.numeric(middayDat[1,]),origin="2016-12-31")

inFileName <- paste(siteName,"_Midday2_varBurn.RData",sep="")
load(inFileName)
var.mat <- as.matrix(var.Burn)
colnames(var.mat)

Tran <- var.mat[,1]
b <- var.mat[,2]
c <- var.mat[,3]
d <- var.mat[,4]
k <- var.mat[,5]
prec <- var.mat[,6]
r <- var.mat[,7]



ci.Midday <- createCI(PFT="SH",var.mat=var.mat,xseq=xseq,doRescale = FALSE)

plot(midday.dates,middayDat[2,],pch=20,ylab="Midday Average",cex.axis=2,cex.lab=2,xlab="")

#polygon(x=c(as.Date("2017-11-30"),as.Date("2017-12-11"),as.Date("2017-12-11"),as.Date("2017-11-30")),y=c(-0.2,-0.2,1.2,1.2),col="#fff7bc",border=NA)
#ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Overall[1,],ci.Overall[3,],col="#d95f02")
ciEnvelope(as.Date(xseq,origin="2016-12-31"),ci.Midday[1,],ci.Midday[3,],col="#1b9e77")

points(midday.dates,middayDat[2,],pch=20)
for(i in 1:ncol(middayDat)){
  print(middayDat[3,i])
  ys <- seq(middayDat[4,i],middayDat[5,i],0.0001)
  #ys <- seq((dat2[2,i]-0.5*dat2[3,i]),(dat2[2,i]+0.5*dat2[3,i]),0.0001)
  #ys <- seq((dat2[2,i]-1.96*dat2[3,i]),(dat2[2,i]+1.96*dat2[3,i]),0.0001)
  lines(rep(midday.dates[i],length(ys)),ys,col="black")
}













##########################
diurnal.files <- dir(pattern=paste("GOES_NDVI_Diurnal",siteName,"_",sep=""))
days <- numeric()
for(i in 1:length(diurnal.files)){
  start <- as.numeric(strsplit(diurnal.files[i],"_")[[1]][4])
  end <- as.numeric(strsplit(diurnal.files[i],"_")[[1]][5])
  dys <- seq(start,end,1)
  days <- c(days,dys)
}
#sort(days)
all.days <- c(seq(1,321,1),seq(347,365,1))
missingDays <- numeric()
for(j in 1:length(all.days)){
  if((!all.days[j] %in% days)){
    missingDays <- c(missingDays,all.days[j])
  }
}
missingDays
length(missingDays)





fileName <- paste("GOES_NDVI_",siteName,"_",startDate,"_",endDate,"_Avg.csv",sep="")



library("ecoforecastR")
library("rjags")


diurnalExp <- function(a,c,k,xseq){
  k <- round(k,digits=1)
  #print(k)
  bk <- which(round(xseq,digits=1)==k)
  #print(bk)
  left <- -a*exp(-1*(xseq[1:bk]-k))+c
  right.xseq <- xseq[(bk+1):length(xseq)]
  right <- -a*exp((right.xseq-k))+c
  #print(length(c(left,right)))
  return(c(left,right))
}


siteName <- "russellSage"
xseq <- seq(0,25,0.1)

outputFileName <- paste(siteName,"_diurnalExamples.pdf",sep="")
pdf(file=outputFileName,width=10,height=5)
par(mfrow=c(2,3))
day.seq <- c("186","201","245","251","278","030")
dy="030"
for(dy in day.seq){
  fileName <- dir(path="dailyNDVI_GOES",pattern=paste(dy,"_varBurn.RData",sep=""))
  print(fileName)
  load(paste("dailyNDVI_GOES/",fileName,sep=""))
  if(as.numeric(dy)<182){
    yr <- "2018"
  }
  else{
    yr <- "2017"
  }
  dat <- read.csv(paste("dailyNDVI_GOES/GOES_Diurnal_",siteName,"_",yr,dy,".csv",sep=""),header=FALSE)
  out.mat <- as.matrix(var.burn)
  a <- out.mat[,1]
  c <- out.mat[,2]
  k <- out.mat[,3]
  ycred <- matrix(0,nrow=10000,ncol=length(xseq))
  for(g in 1:10000){
    Ey <- diurnalExp(a=a[g],c=c[g],k=k[g],xseq=xseq)
    ycred[g,] <- Ey
  }
  ci <- apply(ycred,2,quantile,c(0.025,0.5, 0.975), na.rm= TRUE)
  date <- as.Date(as.numeric(dy),origin=as.Date(paste(as.character(as.numeric(yr)-1),"-12-31",sep="")))
  plot(x=list(),y=list(),main=date,ylim=c(0,1),xlim=c(0,25),ylab="NDVI",xlab="Hour",cex=2.5)
  ciEnvelope(xseq,ci[1,],ci[3,],col="lightBlue")
  lines(xseq,ci[2,],col="black")
  points(as.numeric(dat[3,]),as.numeric(dat[2,]),pch=".")
  abline(v=12,col="red")
}
dev.off()





library("suncalc")
library("lubridate")

getSunlightTimes(date=)
lat <- 32.457
long <- -91.9743
suntimes <- getSunlightTimes(date=as.Date("2017-07-04"),lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),tz = "America/Chicago")
#suntimes <- getSunlightTimes(date=as.Date("2017-07-04"),lat=lat,lon=long,tz = "America/Chicago")
tim <- suntimes$sunrise
dawnTime <- hour(suntimes$nauticalDawn)+(minute(suntimes$nauticalDawn)/60)
duskTime <- hour(suntimes$nauticalDusk)+(minute(suntimes$nauticalDusk)/60)


library("PhenologyBayesModeling")
library("rjags")
library("runjags")
library("MODISTools")
library("ncdf4")
library(plyr)

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
iseq <- c(1)
i <- 2
for(i in iseq){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  URL <- as.character(siteData$URL[i])
  PFT <- as.character(siteData$PFT[i])
  lat <- as.character(siteData$Lat[i])
  long <- as.character(siteData$Long[i])
  TZ <- as.character(siteData$TZ[i])
  #fileName1 <- "GOES_NDVI_HarvardForest_Avg1.csv"
  #fileName2 <- "GOES_NDVI_HarvardForest_Avg2.csv"
  #fileName3 <- "GOES_NDVI_HarvardForest_Avg3.csv"
  #fileName4 <- "GOES_NDVI_HarvardForest_Avg4.csv"
  #fileName5 <- "GOES_NDVI_HarvardForest_Avg5.csv"
  #fileName <- "GOES_NDVI_HarvardForest_Avg.csv"
  #createNDVI_GOES_Avg4(lat=lat,long=long,startDay=152,endDay=546,fileName=fileName4,TZ=TZ)
  
  file1 <- paste("GOES_NDVI_",siteName,"_Avg1.csv",sep="")  
  file2 <- paste("GOES_NDVI_",siteName,"_Avg2.csv",sep="")
  file3 <- paste("GOES_NDVI_",siteName,"_Avg3.csv",sep="")
  file4 <- paste("GOES_NDVI_",siteName,"_Avg4.csv",sep="")
  file5 <- paste("GOES_NDVI_",siteName,"_Avg5.csv",sep="")
  fileALL <- paste("GOES_NDVI_",siteName,"_Avg.csv",sep="")
  dat1 <- read.csv(file1,header=FALSE)
  dat2 <- read.csv(file2,header=FALSE)
  dat3 <- read.csv(file3,header=FALSE)
  dat4 <- read.csv(file4,header=FALSE)
  dat5 <- read.csv(file5,header=FALSE)
  output <- cbind(dat1,dat2,dat3,dat4,dat5)
  write.table(file = fileALL,x=output,col.names=FALSE,row.names=FALSE,sep=",")
  
  
  startDay <- 152
  endDay <- 546
  xseq <- seq(startDay,endDay,1)
  DB.vars <- c("TranF","bF","TranS","bS","c","d","prec","k")
  j.model <- createBayesModel.DB_Avg(dataSource="GOES.NDVI",siteName=siteName,URL="",niter=100000,startDay=152,endDay=546,lat=lat,long=long,TZ=TZ)
  GOES.Avg.md.out <- runMCMC_Model(j.model=j.model,variableNames = DB.vars)
  save(GOES.Avg.md.out,file="HarvardForest_GOES_Avg_varBurn.RData")
  
}

summary(GOES.Avg.md.out)

dat <- read.csv("crappyGOES_NDVI/GOES_NDVI_burns_2017-04-20_2018-03-31_noon.csv",header=FALSE)
plot(as.numeric(dat[1,]),as.numeric(dat[2,]),ylim=c(0,0.5))
load("burns_GOES_varBurn.RData")
out.mat <- as.matrix(GOES.md.out)
colnames(out.mat)

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
iseq <- c(1)
for(i in iseq){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  URL <- as.character(siteData$URL[i])
  PFT <- as.character(siteData$PFT[i])
  lat <- as.character(siteData$Lat[i])
  long <- as.character(siteData$Long[i])
  TZ <- as.character(siteData$TZ[i])
  fileName1 <- "GOES_NDVI_HarvardForest_Avg1.csv"
  fileName2 <- "GOES_NDVI_HarvardForest_Avg2.csv"
  fileName3 <- "GOES_NDVI_HarvardForest_Avg3.csv"
  fileName4 <- "GOES_NDVI_HarvardForest_Avg4.csv"
  fileName5 <- "GOES_NDVI_HarvardForest_Avg5.csv"
  createNDVI_GOES_Avg1(lat=lat,long=long,startDay=startDay,endDay=endDay,fileName=fileName1,TZ=TZ)
  #createNDVI_GOES_Avg1(lat=lat,long=long,startDay=startDay,endDay=endDay,fileName=fileName1,TZ=TZ)
}

GOES <- read.csv("GOES_NDVI_HubbardBrook_Avg.csv",header=FALSE)
GOES2 <- GOES[,colSums(is.na(GOES)) == 0]
#############
GOES.Avg.md.out
load("HarvardForest_GOES_Avg_varBurn.RData")
var.mat<-as.matrix(GOES.Avg.md.out)
GOES.c <- var.mat[,3]
GOES.d <- var.mat[,4]
if(PFT=="DB"){
  GOES.Tran1<-var.mat[,1]
  GOES.Tran2<-var.mat[,2]
  GOES.c <- var.mat[,5]
  GOES.d <- var.mat[,6]
  ci.GOES <- createCI(PFT="DB",var.mat=var.mat,xseq=xseq)
  
  par(mfrow=c(1,1))
  plot(x=list(),y=list(),xlim=c(100,550),ylim=c(-0.2,1.2),ylab="Value",xlab="Day of Year",main=paste(siteName,"GOES"),cex.axis=2,cex.lab=2,cex.main=2)
  lines(xseq,ci.GOES[2,],col="black",lwd=2)
  lines(xseq,ci.GOES[1,],col="black", lty = 2,lwd=2)
  lines(xseq,ci.GOES[3,],col="black", lty = 2,lwd=2)
  abline(v=mean(GOES.Tran1),col="black")
  abline(v=mean(GOES.Tran2),col="black")
  GOES_Days <- as.numeric(GOES2[1,])
  bk <- which(GOES_Days==333)
  print(bk)
  for(i in bk:length(GOES_Days)){
    GOES_Days[i] <- GOES_Days[i]+365
  }
  points(GOES_Days,rescale(c=mean(GOES.c),d=mean(GOES.d),yseq=GOES2[2,]),col="Black",pch=20)
  
}

install.packages("MODISTools")
library("MODISTools")


MODISSubsets(data.frame(lat=Lat,long=Long,start.date=2017,end.date=2018),
             Product="MOD13Q1",Bands="250m_16_days_NDVI",Size=c(1,1),StartDate=TRUE) 

siteName <- "missouriozarks"
startDate <- "2017-04-20"
startDay <- 110
endDay <- 400
startDate <- as.Date(endDay,origin="2017-01-01")
endDate <- "2018-06-27"
lat <- 38.7441
long <- -92.2

directory <- getwd()
metric<- "NDVI"
fileName <- paste(siteName,"_",metric,"_MOD13Q1_",startDate,"_",endDate,".csv",sep="")

if(!file.exists(fileName)){
  mt_subset(product = "MOD13Q1",lat=lat,lon=long,band=paste("250m_16_days_",metric,sep=""),start=startDate,end=endDate,site_name = paste(siteName,"_",metric,sep=""),out_dir = directory,internal=FALSE)
}
dat <- read.csv(fileName,header=TRUE,skip=15)
dat$proc_date[1]
startDateSplit <- c(month=4,day=20,year=2017)
julian(origin=startDate)

x <- numeric()
y <- numeric()

for(i in 1:nrow(dat)){
  y <- c(y,dat$data[i]/10000)
  tmp <- as.Date(dat$calendar_date[i])
  print(tmp)
  x.val <- as.numeric(format(tmp, "%j"))
  print(x.val)
  if(substr(tmp,1,4)=="2018"){
    x.val <- x.val + 365
  }
  x <- c(x,x.val)
}

tmp <- as.Date(dat$calendar_date[22])
format(tmp, "%j",origin=startDateSplit)

MODIS_data(siteName = siteName,lat=lat,long=long,startDay=startDay,endDay=endDay,metric="EVI")


data.MODIS = MODIS_data(siteName=siteName,startDay = startDay,endDay=endDay,lat=lat,long=long,metric=metric)
metric="EVI"


inputFileName <- paste(siteName,"_MODIS_NDVI_varBurn.RData",sep="")
load(inputFileName)

var.mat<-as.matrix(MODIS.N.md.out)
var.mat<-as.matrix(MODIS.E.md.out)


if(PFT=="DB"){
  MODIS.c <- var.mat[,5]
  MODIS.d <- var.mat[,6]
  #MODIS.transF <- mean(-(MODIS.aF/MODIS.bF))
  #MODIS.transS <- mean(-(MODIS.aS/MODIS.bS))
}
else if(PFT=="SH"){
  MODIS.c <- var.mat[,3]
  MODIS.d <- var.mat[,4]
  MODIS.trans <- mean(var.mat[,5])
}
else if(PFT=="EN"){
  MODIS.trans <- mean(var.mat[,2])
}
ci.MODIS <- createCI(PFT=PFT,var.mat=var.mat,xseq=xseq)
print("MODIS Done")

plot(x=list(),y=list(),xlim=c(startDay,endDay),ylim=c(-0.2,1.2),ylab="Value",xlab="Day of Year",main=paste(siteName,"MODIS"),cex.axis=2,cex.lab=2,cex.main=2)
lines(xseq,ci.MODIS[2,],col="red",lwd=2)
lines(xseq,ci.MODIS[1,],col="red", lty = 2,lwd=2)
lines(xseq,ci.MODIS[3,],col="red", lty = 2,lwd=2)
#abline(v=mean(MODIS.trans),col="red")
points(data.MODIS$x,rescale(c=mean(MODIS.c),d=mean(MODIS.d),yseq=data.MODIS$y),col="red",pch=20)

PC.data <- PC_data(siteName=siteName,URL=URL,startDay=startDay,endDay=endDay)






library("PhenologyBayesModeling")
library("rjags")
library("runjags")
library("MODISTools")


graphMCMC_Outputs_withData(outputFileName = "BayesFits_ALL4.pdf",siteFileName = "GOES_Paper_Sites.csv",iseq=c(seq(1,20)),startDaySH=110,endDaySH=455,startDayDB=152,endDayDB=546)
data$x
startDay <- 110
endDay <- 455

#startDay <- 110
#endDay <- 455

siteData <- read.csv("GOES_Paper_Sites.csv",header=TRUE)
i=13

for(i in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  URL <- as.character(siteData$URL[i])
  PFT <- as.character(siteData$PFT[i])
  lat <- as.character(siteData$Lat[i])
  long <- as.character(siteData$Long[i])
  if(PFT=="SH"){
    startDay <- 110
    endDay <- 455
  }
  else if(PFT=="DB"){
    startDay <- 152
    endDay <- 546
  }
  PC_data(siteName=siteName,URL=URL,startDay=startDay,endDay=endDay)
}

data.PC <- PC_data(siteName=siteName,URL=URL,startDay=startDay,endDay=endDay)
data.GOES <- GOES_data(siteName=siteName,lat = lat,long=long,startDay = startDay,endDay = endDay,TZ = TZ)

plot(data.GOES$x,data.GOES$y)
lines(xseq,shrublandYvals(Tran=215,b=-0.13,c=0.25,d=0.15,k=240,r=-0.02,xseq=xseq),col="red")

load("luckyHills_PC_varBurnINITS.RData")
PC.out <- as.matrix(PC.md.out)
PC.a <- mean(PC.out[,1])
PC.b <- mean(PC.out[,2])
PC.c <- mean(PC.out[,3])
PC.d <- mean(PC.out[,4])
PC.k <- mean(PC.out[,5])
PC.r <- mean(PC.out[,7])
#shrublandYvals <- function(Tran,b,c,d,k,r,xseq){

-PC.a/PC.b
data.MODIS = MODIS_data(siteName=siteName,startDay = startDay,endDay=endDay,lat=lat,long=long,metric=metric)
plot(data.MODIS$x,data.MODIS$y,ylim=c(0,1),xlim=c(100,600))
abline(v=365,col="red")
abline(v=550,col="green")
lines(xseq2,pheno.logistic(a=34,b=-0.07,c=0.4,d=0.55,xseq2))
b=
  lines(xseq2,pheno.logistic(a=10,b=-0.04,c=0.4,d=0.55,xseq2))


lines(xseq,deciduousYvals2(aS=34,bS=-0.07,aF=-30,bF=0.11,c=0.4,d=0.55,k=365,xseq=xseq))
xseq1 <- seq(startDay,365,1)
xseq2 <- seq(365,endDay,1)
lines(xseq,deciduousYvals2(aS=30,bS=-0.11,aF=-30,bF=0.11,c=0.4,d=0.55,k=365,xseq=xseq))
plot(data.MODIS$x[14:length(data.MODIS$x)],data.MODIS$y[14:length(data.MODIS$x)],ylim=c(0,1))
lines(xseq2,pheno.logistic(a=34,b=-0.07,c=0.4,d=0.55,xseq2))

##want -a/b =480
a=-480*b
load("Bartlett_MODIS_NDVI_varBurn.RData")
summary(MODIS.N.md.out)

#lines(xseq,deciduousYvals2(aS=aS,bS=bS[g],aF=aF[g],bF=bF[g],c=c[g],d=d[g],k=k[g],xseq=xseq))

deciduousYvals2 <- function(aF,bF,aS,bS,c,d,k,r,xseq){
  bk <- which(xseq==round(k,digits=0))
  print(bk)
  print(xseq[1:bk])
  greendown <- pheno.logistic(a=aF,b=bF,c=c,d=d,xseq[1:bk])
  gu.xseq <- xseq[(bk+1):length(xseq)]
  print(gu.xseq)
  greenup <- pheno.logistic(a=aS,b=bS,c=c,d=d,gu.xseq)
  return(c(greendown,greenup))
}

pheno.logistic <- function(a,b,c,d,xseq){
  return(c/(1 + exp(a+b*xseq))+d)
}







for(i in test){
  print(i)
}



library("ncdf4")
ACM.file <- nc_open("OR_ABI-L1b-RadC-M3C02_G16_s20180010057197_e20180010059569_c20180010100005.nc")
ACM.file <- nc_open("OR_ABI-L1b-RadC-M3C02_G16_s20181411702236_e20181411705009_c20181411705050.nc")
fileCO3 <- nc_open("OR_ABI-L1b-RadC-M3C03_G16_s20180011657198_e20180011659571_c20180011700014.nc")
fileCO2 <- nc_open("OR_ABI-L1b-RadC-M3C02_G16_s20171821657189_e20171821659562_c20171821659598.nc")
structure(fileCO2)


ACM.ind <- getDataIndex(getABI_Index(lat.rd,long.rd),"ACM")
lat <- 42.5378*2*3.14159/360
long <- -72.1715*2*3.14159/360
lat.rd <- as.numeric(lat)*2*pi/360
long.rd <- as.numeric(long)*2*pi/360
Ind2 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion="NEW"),2,orbitVersion="NEW")
Ind3 <- getDataIndex(getABI_Index(lat.rd,long.rd,orbitVersion="NEW"),3,orbitVersion="NEW")

getSpecificNDVI(ind2 = Ind2,ind3 = Ind3,day.time="20180011657")

ACM.ind <- c(-136.7766,415.3229)
ACM.test <- ncvar_get(ACM.file,"BCM")
[ACM.ind[1],ACM.ind[2]]
testgarb <- load("HarvardForest_GOES_varBurn.RData")
summary(testgarb)

library("MODISTools")

lat <- 45.805986
Long <- -90.0791


siteData <- read.csv("GE509_Project_Sites.csv",header=FALSE)

for(i in 1:nrow(siteData)){
  print(i)
  Lat <- as.character(siteData[i,2])
  Long <-as.character(siteData[i,3])
  MODISSubsets(data.frame(lat=Lat,long=Long,start.date=2017,end.date=2018-01-31),
               Product="MOD13Q1",Bands="250m_16_days_NDVI",Size=c(1,1),StartDate=TRUE) 
}


MODISSubsets(data.frame(lat=Lat,long=Long,start.date=2017,end.date=2018),
             Product="MOD13Q1",Bands="250m_16_days_NDVI",Size=c(1,1),StartDate=TRUE) 

logistic <- function(c,a,b,d,t){
  #print(t)
  den <- 1+exp(a+b*t)
  print(den)
  return (c/den + d)
}

exponential <- function(a,r,t){
  return(a*exp(r*t))
}
#EVI
times <- seq(185,365)
c <-0.3 #0.3
d <- 0.25 #0.25
b <- 0.11 #0.11 in fall; -0.11 in spring
a <- -35 #-30 in fall; 30 in spring
vals <- logistic(c,a,b,d,t=times)
plot(times,vals)
hist(rbeta(1000,2,8))

a2 <- 0.3
r <- -0.013
times2 <- seq(0,140)
vals2 <- exponential(a2,r,times2)
plot(times2,vals2)
hist(rgamma(100,5))
```