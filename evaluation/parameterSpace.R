source("simBatchRun.R")

#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- change parameters ----

parlist$description <- NULL
parlist$m <- parlist$n <- 64
parlist$pa <- seq(200, 1200, by=20)
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
parlist$K0 <- 0.2546103
parlist$Kmax <- 2.0461027
parlist$kf <- 1.0
parlist$kc <- 0.6
parlist$useRandomSeed = "T"

#---- calculate parameter space----

Dimnames <- parlist
Dims <- mapply(length, Dimnames)
totalDims <- prod(Dims)
parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)

#---- run Simulation ----
simSpace <- simBatchRun(parameterSpace, folder="/media/Data/eco-hyd/")

#names(simSpace) <- paste0("sim", 1:length(parameterSpace))
#str(simSpace, max.level=1)

#simSpace$sim20$parameter$title

#---- display simulation results ----

# simNo=2
# pos <- which(parameterSpace==simNo, arr.ind=T)
# pos<- as.list(pos)
# 
# variation <- which(names(Dimnames)=="pa")
# pos[[variation]] <- 1:length(Dimnames$pa)
# 
# length(parameterSpace)
# 
# pos <- list(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# parameterSpace[1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1:3 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# mapply(`[`, parameterSpace, pos)
# cat(pos)
# 
# dimnames(parameterSpace)["pa"]

library(R.utils)
#slice.index(parameterSpace, 11)


print(selectedSims <- extract(parameterSpace, run="T",  drop=T))


#eval(parse(text=paste0("sim", selectedSims)))$parameter$pa

result <- selectedSims
for(i in selectedSims){
  result[[i]] <- simSpace[[i]]$postp$coverRatioMedian
}
result

library("fields")
if(is.vector(result)){
  plot(names(result), result, type="l")
}
if(is.matrix(result)){
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

