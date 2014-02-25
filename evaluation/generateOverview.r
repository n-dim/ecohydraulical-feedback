generateOverview <- function(path){
  library(xtable)
  library(foreach)
  library(iterators)
  library(parallel)
  
  source("plotVegetationGrid.r")
  source("clustersize.R")
  source("coverRatio.r")
  grDevices::pdf.options(useDingbats = FALSE);  
  require(knitr)
  opts_chunk$set(fig.path='figure/')
  opts_knit$set(progress = F, verbose = F)
  
  files <- dir(path, pattern='_grids.RData')
  
  oldwd <- getwd()
  
  for(file in files){
    name <- strsplit(file, "_")[[1]][1]
    name <- sub("\\.", "_", name)
    name <- sub(" ", "_", name)
    
    out.folder <- file.path(path, paste0(name, '_overview'))
    dir.create(out.folder, recursive=T)

    out.name <-  paste0(name, '.Rnw')
    file.copy('overview.Rnw', file.path(out.folder,out.name), overwrite=T)
    try({
      setwd(out.folder)
      knit2pdf(out.name, encoding='UTF-8')
    })
    setwd(oldwd)

  }
}
