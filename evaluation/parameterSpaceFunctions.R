library("fields")
library(R.utils)
source("simBatchRun.R")
library("foreach")

plotGridMatrix <- function (sims, title) {
  parold <- par(no.readonly = T)
  try({
  par(mar= c(0,0,0,0), oma=c(4,4,3,0))
  sims <- as.matrix(sims)
  rows <- nrow(sims)
  cols <- ncol(sims)
  layout(t(matrix(1:(rows*cols), nrow=rows))[cols:1,])
  for(i in 1:(rows*cols)){
    sim = sims[i]
    result <- replaceWithResult(sim, simFolder)
    if(all(is.na(result))){
      plot(0,0, type="n", axes=F, asp=1)
      text(0,0, "NA")
    }else{
      image(result$vegetation[[100]], col=gray(9:0/9), breaks=0:10, axes=F, asp=1)
    } 
    #title(main=sims, cex.main=0.8)
    box()
    if(any(i==1:rows)) mtext(rownames(sims)[i], 1, line=1)
    if(i%%(rows)==1) mtext(colnames(sims)[(i+rows-1)/(rows)], 2, line=1)
    #invisible(readline("Enter"))
  }
  mtext(names(dimnames(sims))[1], 1, outer=T, 2)
  mtext(names(dimnames(sims))[2], 2, outer=T, 2)
  mtext(paste(title, collapse=", "), 3, outer=T, 1)
  })
  par(parold)
}

calcParameterSpace <- function (parlist) {
  Dimnames <- parlist
  Dims <- mapply(length, Dimnames)
  totalDims <- prod(Dims)
  message("total Dimensions are = ", totalDims)
  parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)
  return(parameterSpace)
}

createInputFiles <- function (parameterSpace, simFolder) {
  for(simNo in 1:length(parameterSpace)) {
    file <- file.path(simFolder, paste0("input-sim_",simNo,".txt"))
    
    # write input file
    pos <- which(parameterSpace==simNo, arr.ind=T)
    parameterSet <- mapply(`[`, dimnames(parameterSpace), pos)
    write(paste("title =",simNo), file=file)
    for(i in 1:length(parameterSet)){
      par <- recalcParameter(i, parameterSet)
      write(paste(par$name, "=", par$value), file=file, append=T)
    }
  }
}

recalcParameter <- function(index, parameterSet){
  # for input parameters that are not intended in the ecohyd model this function does a transformation
  parname <- names(parameterSet[index])
  parvalue <- parameterSet[index]
  with(as.list(parameterSet), {
    switch(parname,
         Kincrease = {parname <- "Kmax"; parvalue <- as.numeric(K0) + as.numeric(parvalue)},
         KincFrac = {parname <- "Kmax"; parvalue <- as.numeric(K0) * as.numeric(parvalue)}
    )
    return(list(name=parname, value=parvalue))
  })
}

replaceWithResult <- function(simName, simFolder){
  # part can also be "parameter" or "postprocessing"
  grids <- NA
  suppressWarnings(error <- tryCatch(load( paste0(simFolder, "/", simName, "/", simName, "_grids.RData")), error=function(e) NULL, silent=T))
  return(grids)
}



viewParameterSpace <- function(parameterSpace, simFolder, outputParameter, selectiveParameter=list(run=T), plot=T, applyFunction="invisible", randomAverage=T){
  
  print(selectedSims <- extract(parameterSpace, indices=selectiveParameter,  drop=T))
  if(!plot) return(selectedSims)
  # replace simulation number by outputParameter
  result <- selectedSims
  for(i in 1:length(selectedSims)){
    result[[i]] <- NA
    suppressWarnings(error <- tryCatch(load( paste0(simFolder, "/", selectedSims[i], "/", selectedSims[i], "_postprocessing.RData")), error=function(e) NULL, silent=T))
    if(!is.null(error)){
      result[[i]] <- apply(as.matrix(unlist(postprocessing[outputParameter], use.names=F)), 2, FUN=applyFunction)
      rm(postprocessing)  
    }
  }
  if(all(is.na(result))){
    stop("no results yet", call.=F)
  }

  # make an average if parameterSpace uses random seed variation
  if(randomAverage){
    av <- which(names(dimnames(result))=="useRandomSeed")
    if(length(av)!=0){
      dims <- (1:length(dim(result)))[-av]
      resultIQR <- apply(result, dims, function(x) IQR(x, na.rm=T))    
      result <- apply(result, dims, function(x) mean(x, na.rm=T))   
    }
    
  }
  
  print(result)
  
  #--- display results ---
  
  # display in line diagram (if 2d)
  if(is.vector(result) & any(!is.na(result))){
    x <- suppressWarnings(as.numeric(names(result)))
    if(any(is.na(x))) x <- 1:length(x)
    plot(x, result, type="l", main=paste0(names(selectiveParameter), " = ", selectiveParameter), ylab=outputParameter)
    if(randomAverage){
      par(new=T)
      plot(x, resultIQR, col="green", type="l", axes=F, ylab="", xlab="", lty=2)
      axis(4, pretty(resultIQR), col="green")
      mtext("IQR", 4, 3)
    }
  }
  
  # display in shaded picture (if 3d)
  if(is.matrix(result) & any(!is.na(result))){
    x <- suppressWarnings(as.numeric(rownames(result)))
    if(any(is.na(x))) x <- 1:length(x)
    y <- suppressWarnings(as.numeric(colnames(result)))
    if(any(is.na(y))) y <- 1:length(y)
    image.plot(x, y, result, xlab=names(dimnames(result))[1], ylab=names(dimnames(result))[2], main=outputParameter, col=rev(heat.colors(100)), sub="", graphics.reset=F)
    
    # write simulation numbers into grid
   
      text(as.numeric(rep(x, ncol(selectedSims))), as.numeric(rep(y, each=nrow(selectedSims))), selectedSims, cex=.6)

    
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
  invisible(return(selectedSims))
}

printPDFs <- function (withGrids=T, withRandomAverage=F) {
  if(withGrids){
    #--- print grid matrix ---
    asp= ncol(as.matrix(sims))/nrow(as.matrix(sims))
    pdf(paste0(simFolder, "/", names(selectiveParameter), " = ", selectiveParameter, "_gridMatrix.pdf")[1], height=14+7*strheight("x", "inches"), width=14*asp+4*strheight("x", "inches"), onefile=T)
    try({
      plotGridMatrix(sims,title=paste(names(selectiveParameter), "=", selectiveParameter))
    })
    dev.off()
  }
  
  
  #--- print other parameters ----
  
  outputParameters <- c("medianTotalET", "medianTotalBareEvap", "medianTotalDischarge", "medianTotalStore", "medianTotalOutflow", "medianVegDensity", "coverRatioMedian", "wavelength", "wavelength2", "phaseshift", "angularEntropy", "radialEntropy", "Entropy2D", "orientation", "meanWavelength", "meanWavelength2")
  # "wavespeed",
  pdf(paste0(simFolder, "/", names(selectiveParameter), " = ", selectiveParameter, "_parameterPlot.pdf"), onefile=T)
  try({
    for(i in outputParameters){
      viewParameterSpace(parameterSpace, simFolder, outputParameter=i, selectiveParameter=selectiveParameter, plot=T, applyFunction="median", randomAverage=F)
    }
  })
  dev.off()
  
  #--- print with random average ---
  av <- which(names(dimnames(result))=="useRandomSeed")
  
  if(withRandomAverage & length(av)!=0){
    pdf(paste0(simFolder, "/", names(selectiveParameter), " = ", selectiveParameter, "_parameterPlot_randomaverage.pdf"), onefile=T)
    marold <- par()$mar
    par(mar=c(5,4, 4, 4))
    try({
      for(i in outputParameters){
        viewParameterSpace(parameterSpace, simFolder, outputParameter=i, selectiveParameter=selectiveParameter, plot=T, applyFunction="median", randomAverage=T)
      }
    })
    par(mar=marold)
    dev.off()  
  }
}
