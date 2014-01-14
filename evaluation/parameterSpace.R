# for parallel computing:
cores <- 12
library(foreach)
library(doParallel)
registerDoParallel(cores)
source("coverRatio.r")


#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- change parameters ----

parlist$description <- NULL
parlist$m <- parlist$n <- 70
parlist$roughness <- 10^c(-10:13)
parlist$nSteps <- 3

#---- calculate parameter space----

Dimnames <- parlist
#Dimnames <- list(n=c(10,20,30),m=c(10,20,30),A=c("0.1","0.2","0.3"), roughness=c(0.001,0.01,0.1))
Dims <- mapply(length, Dimnames)
#prod(Dims)

parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)
#dimnames(parameterSpace)

#---- create simulation folder ----

folder <- "/media/Data/eco-hyd/"
simFolder <- file.path(folder, paste0("simRun_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(simFolder)

#---- run simulations ----

Data <- new.env()


simSpace <- foreach(simNo=1:length(parameterSpace)) %dopar% {
  # create folder
  simNoFolder <- file.path(simFolder, simNo)
  dir.create(simNoFolder)
  file <- file.path(simNoFolder, paste0("input-sim_",simNo,".txt"))
  
  # write input file
  pos <- which(parameterSpace==simNo, arr.ind=T)
  parameterSet <- mapply(`[`, dimnames(parameterSpace), pos)
  write(paste("title =",simNo), file=file)
  for(i in 1:length(parameterSet)){
    write(paste(names(parameterSet[i]), "=", parameterSet[i]), file=file, append=T)
  }
  
  # run simulation
  system(paste("../model/ecohydModel.out", file,  simNoFolder) )
  
  # read into R format
  Data <-  readCSV(file.path(simNoFolder, paste0(simNo, "_inputParameter.txt")))
  
}

names(simSpace) <- paste0("sim", 1:length(parameterSpace))
str(simSpace, max.level=1)

simSpace$sim20$parameter$title
