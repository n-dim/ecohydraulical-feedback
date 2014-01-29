source("coverRatio.r")
source("readCSV.r")
suppressPackageStartupMessages({
library("fields")
  })


postprocessing <- function(folder=NA) { 
  if(is.na(folder)) stop("please give a folder")
  file <- list.files(folder, "_inputParameter.txt", full.names=T)[1] 
  if(is.na(file)) stop("no file found to read")
  Data <- readCSV(file)
  outFileName <- Data$parameter$title
  
  # postprocessing:
  postprocessing <- NULL
  postprocessing$coverRatio <- coverRatio(Data)
  
  postprocessing$coverRatioMedian <- median(postprocessing$coverRatio)
  
  # save data:
  parameter <- Data$parameter
  save(parameter, file=file.path(folder, paste0(outFileName, "_parameter.RData")))  
  grids <- Data$rasters
  save(grids, file=file.path(folder, paste0(outFileName, "_grids.RData")))        
  save(postprocessing, file=file.path(folder, paste0(outFileName, "_postprocessing.RData")))
  #save(list=outFileName, file=file.path(folder,paste(outFileName, "_grids.RData", sep="")))  
}

