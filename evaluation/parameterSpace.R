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

#---- calculate parameter space----

Dimnames <- parlist
Dims <- mapply(length, Dimnames)
totalDims <- prod(Dims)
parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)

#---- create simulation folder ----
simFolder <- file.path(folder, paste0("simRun_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(simFolder)


#---- create input files ----
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

#---- run Simulation ----
system(paste("./parallelBatchRun.sh", simFolder), wait=F)


#---- explore parameter space ----
print(selectedSims <- extract(parameterSpace, run="T",  drop=T))

result <- selectedSims

for(i in selectedSims){
  result[[i]] <- NA
  suppressWarnings(error <- tryCatch(load( paste0(simFolder, "/", i, "/", i, "_postprocessing.RData")), error=function(e) NULL, silent=T))
  if(!is.null(error)){
    result[[i]] <- postprocessing$coverRatioMedian
    rm(postprocessing)  
  }
}
result

#--- display results ---

# display in line diagram (if 2d)
if(is.vector(result) & any(!is.na(result))){
  plot(names(result), result, type="l")
}

# display in shaded picture (if 3d)
if(is.matrix(result) & any(!is.na(result))){
  image.plot(as.numeric(rownames(result)), as.numeric(colnames(result)), result, las=1, xlab=names(dimnames(result))[1], ylab=names(dimnames(result))[2], legend.lab="median vegetation cover ratio", col=rev(heat.colors(100)), sub="", graphics.reset=T)

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


for(i in 1:Dimnames$nSteps[1]){
  image.plot(simSpace[[extract(parameterSpace, pa="1200", Emax=1,  drop=T)]]$rasters$vegetation[[i]], col=gray(8:0/8), breaks=0:9-0.5)
  invisible(readline(i))
}


points <- identify(names(result), result, labels=names(result))
points(names(result)[points], result[points])

