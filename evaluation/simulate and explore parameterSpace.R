source("parameterSpace.R")

folder <- "/media/Data/eco-hyd/" #simulation output folder

#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- change parameters ----

parlist$description <- NULL
parlist$m <- parlist$n <- 60
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
#parlist$K0  <- c(0.005, 0.05, 0.5, 1, 2)
parlist$K0 <- seq(0.5, 1, by=0.2)
#parlist$Kmax <- 2.0461027
#parlist$Kmax <- seq(0.5, 6, by=0.5)
parlist$Kmax <- NULL # to calc Kmax as K0 + Kincrease or Kmax * KincFrac
#parlist$Kincrease <- seq(0.1, 2, by=0.2)
parlist$Kincrease <- NULL
#parlist$Kincrease <- c(1, 2, 3, 4)
#parlist$KincFrac <- 1 + 2^c(-2,-1,0,1,2)
parlist$KincFrac <- 1 + 2^c(2,3,4)
parlist$kf <- 1.0
parlist$kc <- 0.6
parlist$useRandomSeed = "T"
#parlist$m <- parlist$n <- 100
#parlist$pa <- seq(200, 1200, by=20)
#parlist$Kmax <- 2.0461027 + seq(0.1, 2, by=0.1)

#--- run no 5----
parlist$KincFrac <- seq(5,15, 2)
parlist$K0 <- 0.5

#--- run no 6 ---
parlist$pa <- c(400,600, 800)
parlist$KincFrac <- c(4,9,15,16,17) 

#--- run no 7 ---
parlist$pa <- 600
parlist$KincFrac <- c(2,3,4,9,15,16,17,18) 
parlist$K0 <- c(0.025, 0.5)

#--- run no 8 ---
parlist$pa <- 600
parlist$KincFrac <- c(2,3,4,9,15,16,17,18) 
parlist$K0 <- seq(0.06, 0.5, by=0.04)

#--- run no 9pre ---
parlist$KincFrac <- c(4,4.5,5,6,9) 
parlist$K0 <- 0.26

#--- run no 9 ---
parlist$KincFrac <- c(4,5,6,9,15,16,17,18,19,20) 
parlist$K0 <- seq(0.02, 0.5, by=0.04)



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
par(mfrow=c(1,1), mar=c(4,4,3,1), oma=rep(0,4))

outputParameters <- c("medianTotalET", "medianTotalBareEvap", "medianTotalDischarge", "medianTotalStore", "medianTotalOutflow", "medianVegDensity", "coverRatioMedian", "wavenumber", "wavelength2", "angularEntropy", "radialEntropy", "Entropy2D", "orientation")

selectiveParameter <- list( K0=-1:-4)
outputParameter <- "medianVegDensity"
outputParameter <- "wavelength2"
sims <- viewParameterSpace(parameterSpace, simFolder, outputParameter=outputParameter, selectiveParameter=selectiveParameter, plot=T, applyFunction="median")

#--- print grid matrix ---
asp= nrow(sims)/ncol(sims)
library(Cairo)
pdf(paste0(simFolder, "/", names(selectiveParameter), " = ", selectiveParameter, "_gridMatrix.pdf")[1], height=14+7*strheight("x", "inches"), width=14*asp+4*strheight("x", "inches"), onefile=T)

plotGridMatrix(sims,title=paste(names(selectiveParameter), "=", selectiveParameter))

dev.off()

#--- print other parameters ----
pdf(paste0(simFolder, "/", names(selectiveParameter), " = ", selectiveParameter, "_parameterPlot.pdf"), onefile=T)

for(i in outputParameters){
  viewParameterSpace(parameterSpace, simFolder, outputParameter=i, selectiveParameter=selectiveParameter, plot=T, applyFunction="median")
}

dev.off()

