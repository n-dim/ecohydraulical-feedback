source("coverRatio.r")
source("readCSV.r")
suppressPackageStartupMessages({
library("fields")
  })


postprocessing <- function(folder=NA) { 
  if(is.na(folder)) stop("please give a folder")
  file <- list.files(folder, "_inputParameter.txt", full.names=T)
  if(length(file)>1) message("only the first output file is processed")
  if(length(file)==0) stop("no file to postprocess")
  file <- file[1]
  if(is.na(file)) stop("no file found to read")
  Data <- readCSV(file)
  outFileName <- Data$parameter$title
  
  # postprocessing:
  postprocessing <- NULL
  postprocessing$coverRatio <- coverRatio(Data)
  postprocessing$coverRatioMedian <- median(postprocessing$coverRatio)
  
  postprocessing$medianTotalET <- median(Data$Summary$totalET)  
  postprocessing$medianVegDensity <- median(Data$Summary$vegDensity)
  postprocessing$medianTotalBE <- median(Data$Summary$totalBE)
  postprocessing$medianTotalStore <- median(Data$Summary$totalStore)
  postprocessing$medianTotalDischarge <- median(Data$Summary$totalDischarge)
  postprocessing$medianTotalOutflow <- median(Data$Summary$totalOutflow)
    
  # save data:
  parameter <- Data$parameter
  save(parameter, file=file.path(folder, paste0(outFileName, "_parameter.RData")))  
  grids <- Data$rasters
  save(grids, file=file.path(folder, paste0(outFileName, "_grids.RData")))        
  save(postprocessing, file=file.path(folder, paste0(outFileName, "_postprocessing.RData")))
  #save(list=outFileName, file=file.path(folder,paste(outFileName, "_grids.RData", sep="")))  
}

