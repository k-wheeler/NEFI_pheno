#days <- c(seq(1,333),seq(348,365))
days <- c(seq(98,120),seq(221-232),seq(360-365))
missingC02_0 <- character()
missingC03_0 <- character()
missingACM_0 <- character()
missingC02_1 <- character()
missingC03_1 <- character()
missingACM_1 <- character()
missingC02_2 <- character()
missingC03_2 <- character()
missingACM_2 <- character()
for(i in 1:length(days)){
  if(as.numeric(days[i]) < 182){
    year <- 2018
  }
  else{
    year <- 2017
  }
  if(as.numeric(days[i])<10){
    days[i] <- as.character(paste("00",as.character(days[i]),sep=""))
  }
  else if(as.numeric(days[i])<100){
    days[i] <- as.character(paste("0",as.character(days[i]),sep=""))
  }
  else{
    days[i] <- as.character(days[i])
  }

  print(paste("s",year,days[i],sep=""))
  ##0:
  fileStr <- paste("s",year,days[i],"0",sep="")
  C02.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C02"))
  C03.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C03"))
  ACM.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C02"))
  
  if(length(C02.files)==0){
    missingC02_0 <- c(missingC02_0,fileStr)
  }
  if(length(C03.files)==0){
    missingC03_0 <- c(missingC03_0,fileStr)
  }
  if(length(ACM.files)==0){
    missingACM_0 <- c(missingACM_0,fileStr)
  }

  ##1:
  fileStr <- paste("s",year,days[i],"1",sep="")
  C02.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C02"))
  C03.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C03"))
  ACM.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C02"))
  
  if(length(C02.files)==0){
    missingC02_1 <- c(missingC02_1,fileStr)
  }
  if(length(C03.files)==0){
    missingC03_1 <- c(missingC03_1,fileStr)
  }
  if(length(ACM.files)==0){
    missingACM_1 <- c(missingACM_1,fileStr)
  }
  
  ##2:
  fileStr <- paste("s",year,days[i],"1",sep="")
  C02.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C02"))
  C03.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C03"))
  ACM.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C02"))
  
  if(length(C02.files)==0){
    missingC02_2 <- c(missingC02_2,fileStr)
  }
  if(length(C03.files)==0){
    missingC03_2 <- c(missingC03_2,fileStr)
  }
  if(length(ACM.files)==0){
    missingACM_2 <- c(missingACM_2,fileStr)
  }
}
print("C02_0:")
print(missingC02_0)
print("C03_0:")
print(missingC03_0)
print("ACM_0:")
print(missingACM_0)
print("C02_1:")
print(missingC02_1)
print("C03_1:")
print(missingC03_1)
print("ACM_1:")
print(missingACM_1)
print("C02_2:")
print(missingC02_2)
print("C03_2:")
print(missingC03_2)
print("ACM_2:")
print(missingACM_2)
