source("coverRatio.r")
source("readCSV.r")
suppressPackageStartupMessages({
library("fields")
  })


postprocessing <- function(outputFolder=NA) { 
  if(is.na(outputFolder)) stop("please give a folder")
  file <- list.files(outputFolder, "_inputParameter.txt", full.names=T)
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
  postprocessing$medianTotalBareEvap <- median(Data$Summary$totalBE)
  postprocessing$medianTotalStore <- median(Data$Summary$totalStore)
  postprocessing$medianTotalDischarge <- median(Data$Summary$totalDischarge)
  postprocessing$medianTotalOutflow <- median(Data$Summary$totalOutflow)
  
  # wavelength and orientation:
  source("pattern recognition/analyseSpectrum.R")
  source("pattern recognition/extract wavelength.R")
  Spectrum <- list(NULL)
  for (t in 1: Data$parameter$nSteps){
    grid <- Data$rasters$vegetation[[t]]
    #grid <- grid[rep(1:nrow(grid), each=1), rep(1:ncol(grid), each=1) ]
    Spectrum <- analyseSpectrum(image=grid, F, F, F, F)
    postprocessing$wavenumber[t] <- Spectrum$radial$mids[which.max(Spectrum$radial$power)]
    postprocessing$wavelength[t] <- nrow(grid)/Spectrum$radial$mids[which.max(Spectrum$radial$power)]
    postprocessing$orientation[t] <- Spectrum$angular$mids[which.max(Spectrum$angular$power)]
    postprocessing$radialEntropy[t] <- Spectrum$entropy_radial
    postprocessing$angularEntropy[t] <- Spectrum$entropy_angular
    postprocessing$Entropy2D[t] <- Spectrum$entropy_2D
    Temp <-  extractWaveLength(bild=Data$rasters$vegetation[[t]], F, F, smoothing=1.5) 
    postprocessing$wavelength2[t] <- Temp $ wavelength
    postprocessing$phaseshift[t] <- Temp $ phaseshift
    rm(Temp)
    rm(grid)
  }
  postprocessing$wavespeed <- c(NA, diff(postprocessing$phaseshift))
  
  # save data:
  parameter <- Data$parameter
  save(parameter, file=file.path(outputFolder, paste0(outFileName, "_parameter.RData")))  
  grids <- Data$rasters
  save(grids, file=file.path(outputFolder, paste0(outFileName, "_grids.RData")))        
  save(postprocessing, file=file.path(outputFolder, paste0(outFileName, "_postprocessing.RData")))
  #save(list=outFileName, file=file.path(outputFolder,paste(outFileName, "_grids.RData", sep="")))  
}

