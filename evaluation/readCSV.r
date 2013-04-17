readCSV <- function(file=NA) {
  # Choose File with input parameters to read in first
  if(is.na(file)) file <- file.choose()
  Temp <- read.table(file, sep="=", colClasses="character", strip.white=T, nrows=35)
  # reform as list
  parameter <- as.list(Temp[,2])
  # give list elements the right names
  names(parameter) <- Temp[,1]
  formats <- c(rep("character", 2), "logical", rep("numeric", 25), rep("logical", 7))
  for(i in 1:length(parameter)){
    parameter[i] <- as(parameter[i], formats[i])
  }
  
  attach(parameter)
  
  #only read rasters if simulation was run
  if(run == T){ 
    # prepare raster read in
    rasterSets <- c("bareE", "discharge", "eTActual", "flowdirections", "flowResistance", "store", "topography", "vegetation")
    rasters <- vector("list", length(rasterSets))
    # read in rasters
    for(i in 1:length(rasterSets)){
      connection <- file(paste(dirname(file),"/", title, "_", rasterSets[i], ".csv", sep=""), open="rt")
      for(j in 1:nSteps){
        Temp <- as.matrix(read.table(file=connection, sep=";", skip=1, nrow=m, header=F))
        colnames(Temp) <- NULL  
        rasters[[i]][[j]] <- Temp[1:n,1:m]
      }
      close(connection)
    }
    # rename it
    names(rasters) <- rasterSets
    
    # read summary file
    Summary <- read.table(file=paste(dirname(file),"/", title, "_", "SummaryResults.csv", sep=""), sep=";", header=T)[,-8] #last column (#8) is empty
    
    message('read input from parameterset "', parameter$title, '"')
    detach('parameter')
    return(list(parameter=parameter, Summary=Summary, rasters=rasters))
  
  }else{
    message('no data to read for parameterset "', parameter$title, '" (no simulation run)')
    detach("parameter")
    return(list(parameter=parameter))
  }
  
  
}