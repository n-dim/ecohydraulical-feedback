source("parameterSpace.R")

folder <- "/media/Data/eco-hyd/" #simulation output folder

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