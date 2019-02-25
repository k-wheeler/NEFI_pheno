#days <- c(seq(1,333),seq(348,365))
days <- c(seq(183,365),seq(1,37))
missingC02 <- character()
missingC03 <- character()
missingC05 <- character()
missingACM <- character()

for(i in 1:length(days)){
  if(as.numeric(days[i]) > 182){
    year <- 2018
  }
  else{
    year <- 2019
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

  fileStr <- paste("s",year,days[i],sep="")
  C02.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C02"))
  C03.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C03"))
  C05.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="M3C05"))
  ACM.files <- intersect(dir(path="GOES_Data2017",pattern=fileStr),dir(path="GOES_Data2017",pattern="ACM"))
  
  if(length(C02.files)==0){
    missingC02 <- c(missingC02,fileStr)
  }
  if(length(C03.files)==0){
    missingC03 <- c(missingC03,fileStr)
  }
  if(length(C05.files)==0){
    missingC05 <- c(missingC05,fileStr)
  }
  if(length(ACM.files)==0){
    missingACM <- c(missingACM,fileStr)
  }
}

print("C02:")
print(missingC02)
print("C03:")
print(missingC03)
print("C05:")
print(missingC05)
print("ACM:")
print(missingACM)

output <- cbind(missingC02,missingC03,missingC05,missingACM)
table.write(file="missingGOESFiles.csv",output,row.names=FALSE,col.names=TRUE)


