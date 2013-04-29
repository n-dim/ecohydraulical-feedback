source("readCSV.r")
source("CSVtoRData.r")

readAllCSV <- function(path){
  files <- dir(path, pattern="inputParameter.txt")
  wdOld <- getwd()
  setwd(path)
  try(
    for (file in files){
      message("converting files associated with ", file, " to .RData")
      CSVtoRData(file)  
    }
  )
  setwd(wdOld)
  
}
