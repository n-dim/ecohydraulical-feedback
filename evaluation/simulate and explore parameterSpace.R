source("parameterSpace.R")

folder <- "/media/Data/eco-hyd/" #simulation output folder

#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- change parameters ----

parlist$description <- NULL
parlist$m <- parlist$n <- 100
#parlist$pa <- seq(200, 1100, by=50)
parlist$pa <- c(400,600, 800,1200)
#parlist$pa <- 600
parlist$roughness <- 0.0000001
parlist$Emax <- 876
#parlist$Emax <- seq(200, 2000, by=100)
#parlist$roughness <- 0
#parlist$nSteps <- 100
parlist$BCs <- "-3\t-3\t0\t0"
parlist$simErosion <- "F"
#parlist$np <- 150
parlist$np <- 0 # np is then calculated as pa/4
parlist$nSteps <- 100
#parlist$K0 <- 0.2546103 + seq(0.0, 10, by=1)
parlist$K0  <- c(0.005, 0.05, 0.5, 1, 2)
#parlist$Kmax <- 2.0461027
#parlist$Kmax <- seq(0.5, 6, by=0.5)
parlist$Kmax <- 0 # to calc Kmax as K0 + Kincrease or Kmax * KincFrac
#parlist$Kincrease <- seq(0.1, 2, by=0.2)
parlist$Kincrease <- 0
#parlist$Kincrease <- c(1, 2, 3, 4)
parlist$KincFrac <- 1 + 2^c(-2,-1,0,1,2)
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
par(mfrow=c(1,1), mar=c(4,4,3,1))

#   "coverRatioMedian"
#   "medianTotalET"
#   "medianTotalDischarge"
#   "medianTotalStore"
#   "medianVegDensity"
#   "medianVegDensity"
#   "medianTotalBareEvap"
#   "medianTotalOutflow"

selectiveParameter <- list(pa="600")
outputParameter <- "medianVegDensity"
sims <- viewParameterSpace(parameterSpace, simFolder, outputParameter=outputParameter, selectiveParameter=selectiveParameter, plot=F)
plotGridMatrix(sims,title=paste(names(selectiveParameter), "=", selectiveParameter))


#------------------
viewParameterSpace(OldParameterSpace, "/media/Data/eco-hyd//simRun_2014-02-14_14-58-08", "coverRatioMedian")


########################################


for(i in 1:Dimnames$nSteps[1]){
  image.plot(simSpace[[extract(parameterSpace, pa="1200", Emax=1,  drop=T)]]$rasters$vegetation[[i]], col=gray(8:0/8), breaks=0:9-0.5)
  invisible(readline(i))
}


points <- identify(names(result), result, labels=names(result))
points(names(result)[points], result[points])

