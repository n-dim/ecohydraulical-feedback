source("parameterSpaceFunctions.R")
#source("simulate_ParameterSpace.R")


parameterSpace <- file.choose()
simFolder <- dirname(parameterSpace)
source(parameterSpace)

system(paste("nautilus", simFolder, "&"))

#---- explore parameter space ----
par(mfrow=c(1,1), mar=c(4,4,3,1), oma=rep(0,4))


selectiveParameter <- list(run=T)
outputParameter <- "medianVegDensity"
outputParameter <- "wavelength2"
outputParameter <- "wavenumber"
outputParameter <- "wavelength"
outputParameter <- "wavespeed"
outputParameter <- "medianTotalET"
outputParameter <- "phaseshift"
outputParameter <- "meanWavelength"

sims <- viewParameterSpace(parameterSpace, simFolder, outputParameter=outputParameter, selectiveParameter=selectiveParameter, plot=T, applyFunction="median", randomAverage=F)

# print pdfs with all overview graphics into simulation Folder:
printPDFs(withRandomAverage=T, withGrids=T)
