library("fields")
library(R.utils)
source("simBatchRun.R")
library("foreach")

folder <- "/media/Data/eco-hyd/" #simulation output folder

#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- change parameters ----

parlist$description <- NULL
parlist$m <- parlist$n <- 20
parlist$pa <- seq(200, 1100, by=50)
#parlist$pa <- 600
parlist$roughness <- 0.0000001
parlist$Emax <- 876
#parlist$Emax <- seq(200, 2000, by=100)
#parlist$roughness <- 0
#parlist$nSteps <- 100
parlist$BCs <- "-3\t-3\t0\t0"
parlist$simErosion <- "F"
parlist$np <- 150
parlist$nSteps <- 100
parlist$K0 <- 0.2546103 + seq(0.1, 0.5, by=0.1)
parlist$Kmax <- 2.0461027
parlist$kf <- 1.0
parlist$kc <- 0.6
parlist$useRandomSeed = "T"


#parlist$m <- parlist$n <- 100
#parlist$pa <- seq(200, 1200, by=20)
#parlist$Kmax <- 2.0461027 + seq(0.1, 2, by=0.1)

#---- create parameter space ----
parameterSpace <- calcParameterSpace(parlist)

#---- create simulation folder ----
simFolder <- file.path(folder, paste0("simRun_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(simFolder)

#---- create input files ----
createInputFiles(parameterSpace, simFolder)
# dump parameter Space
dump(list="parameterSpace", file=file.path(simFolder, "parameterSpace.R"))

#---- run Simulation ----
system(paste("cd", simFolder, "&&", "parallelBatchRun.sh", simFolder), wait=F)


#---- explore parameter space ----
par(mfrow=c(3,1), mar=c(4.5,4,0.1,1.1))

viewParameterSpace(parameterSpace, simFolder, "coverRatioMedian")
viewParameterSpace(parameterSpace, simFolder, "medianTotalET")
viewParameterSpace(parameterSpace, simFolder, "medianTotalDischarge")
viewParameterSpace(parameterSpace, simFolder, "medianTotalStore")
viewParameterSpace(parameterSpace, simFolder, "medianVegDensity")
viewParameterSpace(parameterSpace, simFolder, "medianTotalBE")
#viewParameterSpace(parameterSpace, simFolder, "medianTotalOutflow")

viewParameterSpace(OldParameterSpace, "/media/Data/eco-hyd//simRun_2014-02-14_14-58-08", "coverRatioMedian")

#--- function definitions ----

calcParameterSpace <- function (parlist) {
  #---- calculate parameter space----
  
  Dimnames <- parlist
  Dims <- mapply(length, Dimnames)
  totalDims <- prod(Dims)
  message("total Dimensions are = ", totalDims)
  parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)
  return(parameterSpace)
}

createInputFiles <- function (parameterSpace, simFolder) {
  simSpace <- foreach(simNo=1:length(parameterSpace)) %do% {
    file <- file.path(simFolder, paste0("input-sim_",simNo,".txt"))
    
    # write input file
    pos <- which(parameterSpace==simNo, arr.ind=T)
    parameterSet <- mapply(`[`, dimnames(parameterSpace), pos)
    write(paste("title =",simNo), file=file)
    for(i in 1:length(parameterSet)){
      write(paste(names(parameterSet[i]), "=", parameterSet[i]), file=file, append=T)
    }
  }
}


viewParameterSpace <- function(parameterSpace, simFolder, outputParameter, selectiveParameter=list(run=T)){
  
  print(selectedSims <- extract(parameterSpace, indices=selectiveParameter,  drop=T))
  
  # replace simulation number by outputParameter
  result <- selectedSims
  for(i in selectedSims){
    result[[i]] <- NA
    suppressWarnings(error <- tryCatch(load( paste0(simFolder, "/", i, "/", i, "_postprocessing.RData")), error=function(e) NULL, silent=T))
    if(!is.null(error)){
      result[[i]] <- unlist(postprocessing[outputParameter])
      rm(postprocessing)  
    }
  }
  print(result)
  
  #--- display results ---
  
  # display in line diagram (if 2d)
  if(is.vector(result) & any(!is.na(result))){
    plot(names(result), result, type="l", main=paste0(names(selectiveParameter), " = ", selectiveParameter), ylab=outputParameter)
  }
  
  # display in shaded picture (if 3d)
  if(is.matrix(result) & any(!is.na(result))){
    image.plot(as.numeric(rownames(result)), as.numeric(colnames(result)), result, xlab=names(dimnames(result))[1], ylab=names(dimnames(result))[2], legend.lab=outputParameter, col=rev(heat.colors(100)), sub="", graphics.reset=T)
    
    # other parameters in figure label:
    pos <- which(parameterSpace==selectedSims[1], arr.ind=T)
    annotation <- mapply(`[`, dimnames(parameterSpace), pos)
    for(i in 1:2){
      annotation[which(names(annotation)==names(dimnames(selectedSims))[i])] <- paste0(c("x","y")[i], "-axis")
    }
    paste(names((annotation)), "=", annotation, collapse=", ")
    
  }
  if(length(dim(selectedSims))==3){
    dim(selectedSims)[3]
  }
}




########################################


for(i in 1:Dimnames$nSteps[1]){
  image.plot(simSpace[[extract(parameterSpace, pa="1200", Emax=1,  drop=T)]]$rasters$vegetation[[i]], col=gray(8:0/8), breaks=0:9-0.5)
  invisible(readline(i))
}


points <- identify(names(result), result, labels=names(result))
points(names(result)[points], result[points])

