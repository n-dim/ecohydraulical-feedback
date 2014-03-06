source("parameterSpaceFunctions.R")
source("simulate_ParameterSpace.R")

system(paste("nautilus", simFolder))

#---- explore parameter space ----
par(mfrow=c(1,1), mar=c(4,4,3,1), oma=rep(0,4))

outputParameters <- c("medianTotalET", "medianTotalBareEvap", "medianTotalDischarge", "medianTotalStore", "medianTotalOutflow", "medianVegDensity", "coverRatioMedian", "wavenumber", "wavelength", "wavelength2", "phaseshift", "wavespeed", "angularEntropy", "radialEntropy", "Entropy2D", "orientation")

selectiveParameter <- list(run=T)
outputParameter <- "medianVegDensity"
outputParameter <- "wavelength2"
outputParameter <- "wavenumber"
outputParameter <- "wavelength"
outputParameter <- "wavespeed"
sims <- viewParameterSpace(parameterSpace, simFolder, outputParameter=outputParameter, selectiveParameter=selectiveParameter, plot=T, applyFunction="median", randomAverage=F)

#--- print grid matrix ---
asp= nrow(sims)/ncol(sims)
pdf(paste0(simFolder, "/", names(selectiveParameter), " = ", selectiveParameter, "_gridMatrix.pdf")[1], height=14+7*strheight("x", "inches"), width=14*asp+4*strheight("x", "inches"), onefile=T)

plotGridMatrix(sims,title=paste(names(selectiveParameter), "=", selectiveParameter))

dev.off()

#--- print other parameters ----
pdf(paste0(simFolder, "/", names(selectiveParameter), " = ", selectiveParameter, "_parameterPlot.pdf"), onefile=T)

for(i in outputParameters){
  viewParameterSpace(parameterSpace, simFolder, outputParameter=i, selectiveParameter=selectiveParameter, plot=T, applyFunction="median", randomAverage=F)
}

dev.off()

