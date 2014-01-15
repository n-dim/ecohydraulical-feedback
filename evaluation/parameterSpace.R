# for parallel computing:
cores <- 12
library(foreach)
library(doParallel)
registerDoParallel(cores)
source("coverRatio.r")
source("readCSV.r")


#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- change parameters ----

parlist$description <- NULL
parlist$m <- parlist$n <- 10
parlist$pa <- seq(200, 800, by=20)
parlist$roughness <- c(0)
#parlist$Emax <- 876
parlist$Emax <- seq(300, 2000, length.out=10)
#parlist$roughness <- 0
#parlist$nSteps <- 100

#---- calculate parameter space----

Dimnames <- parlist
Dims <- mapply(length, Dimnames)
totalDims <- prod(Dims)

parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)
#dimnames(parameterSpace)

#---- create simulation folder ----

folder <- "/media/Data/eco-hyd/"
simFolder <- file.path(folder, paste0("simRun_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
dir.create(simFolder)

#---- run simulations ----
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
  system(paste("../model/ecohydModel.out", file,  simNoFolder, ">", paste0(simNoFolder, "/", simNo, "_log.txt")))
  
  # read into R format
  Data <-  readCSV(file.path(simNoFolder, paste0(simNo, "_inputParameter.txt")))
  
  Data$postp$coverRatio <- coverRatio(Data)
  
  Data$postp$coverRatioMedian <- median(Data$postp$coverRatio)
  
  # save Data
  outFileName <- paste0("sim", simNo)
  assign(outFileName, Data)
  save(list=outFileName, file=file.path(simNoFolder,paste(simNo, "_grids.RData", sep="")))  
  
  cat("finished simulation no", simNo, "of", length(parameterSpace), file=file.path(simFolder,"proceeding.txt"), fill=T)
  
  Data
}

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

if(is.vector(result)){
  plot(names(result), result)
}
if(is.matrix(result)){
  image.plot(as.numeric(rownames(result)), as.numeric(colnames(result)), result, las=1, xlab=names(dimnames(result))[1], ylab=names(dimnames(result))[2], legend.lab="median vegetation cover ratio", col=rev(heat.colors(100)), sub=paste(1:10, collapse="")  )
}
if(length(dim(selectedSims))==3){
  dim(selectedSims)[3]
}

pos <- which(parameterSpace==selectedSims[1], arr.ind=T)
annotation <- mapply(`[`, dimnames(parameterSpace), pos)
annotation[which(names(annotation)=="Emax")] <- "x-axis"
paste(names((annotation)), "=", annotation, collapse="; ")



cat(rep(1, length(parameterSpace)))

