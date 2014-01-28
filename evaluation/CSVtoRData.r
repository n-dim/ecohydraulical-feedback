source("readCSV.r")

CSVtoRData <- function(file=NA) {  #choose ..._inputParameter.txt as input File
  if(is.na(file)) file <- file.choose()
  folder <- dirname(file)
  Data <- readCSV(file)
  outFileName <- Data$parameter$title
  assign(outFileName, Data)
  save(list=outFileName, file=file.path(folder,paste(outFileName, "_grids.RData", sep="")))
  #save(list=c(parameter=Data$parameter), file=paste(outFileName), sep="", "_inputParameter.RData")
  
  #dump(ls(Data), file=paste(outFileName, ".RData", sep=""))
}
