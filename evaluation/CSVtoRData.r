source("readCSV.r")

CSVtoRData <- function(file) {  #choose ..._inputParameter.txt as input File
  if(is.na(file)) file <- file.choose()
  Data <- readCSV(file)
  outFileName <- Data$parameter$title
  assign(outFileName, Data)
  rm(Data)
  save(list=outFileName, file=paste(outFileName, ".RData", sep=""))
  #dump(ls(Data), file=paste(outFileName, ".RData", sep=""))
}
