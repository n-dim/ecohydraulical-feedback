library("fields")
library(R.utils)
source("simBatchRun.R")
library("foreach")

plotGridMatrix <- function (sims, title) {
  parold <- par(no.readonly = T)
  try({
  par(mar= rep(0.5,4), oma=c(4,4,3,0))
  rows <- nrow(sims)
  cols <- ncol(sims)
  layout(t(matrix(1:(rows*cols), nrow=rows))[cols:1,])
  for(i in 1:(rows*cols)){
    sim = sims[i]
    result <- replaceWithResult(sim, simFolder)
    if(is.na(result[1])){
      plot(0,0, type="n", axes=F)
      text(0,0, "NA")
      box()
    }else{
      image(result$vegetation[[100]], col=gray(9:0/9), breaks=0:10, axes=F, asp=1, main=sim, cex.main=0.7)
      box()
    } 
    
    if(any(i==1:rows)) mtext(rownames(sims)[i], 1, line=1)
    if(i%%rows==1) mtext(colnames(sims)[(i-1+cols)/cols], 2, line=1)
    #invisible(readline("Enter"))
  }
  mtext(names(dimnames(sims))[1], 1, outer=T, 2)
  mtext(names(dimnames(sims))[2], 2, outer=T, 2)
  mtext(title, 3, outer=T, 1)
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
}

replaceWithResult <- function(simName, simFolder){
  # part can also be "parameter" or "postprocessing"
  grids <- NA
  suppressWarnings(error <- tryCatch(load( paste0(simFolder, "/", simName, "/", simName, "_grids.RData")), error=function(e) NULL, silent=T))
  return(grids)
}



viewParameterSpace <- function(parameterSpace, simFolder, outputParameter, selectiveParameter=list(run=T), plot=T){
  
  print(selectedSims <- extract(parameterSpace, indices=selectiveParameter,  drop=T))
  if(!plot) return(selectedSims)
  # replace simulation number by outputParameter
  result <- selectedSims
  for(i in 1:length(selectedSims)){
    result[[i]] <- NA
    suppressWarnings(error <- tryCatch(load( paste0(simFolder, "/", selectedSims[i], "/", selectedSims[i], "_postprocessing.RData")), error=function(e) NULL, silent=T))
    if(!is.null(error)){
      result[[i]] <- unlist(postprocessing[outputParameter])
      rm(postprocessing)  
    }
  }
  if(all(is.na(result))){
    stop("no results yet", call.=F)
  }
  print(result)
  
  #--- display results ---
  
  # display in line diagram (if 2d)
  if(is.vector(result) & any(!is.na(result))){
    plot(names(result), result, type="l", main=paste0(names(selectiveParameter), " = ", selectiveParameter), ylab=outputParameter)
  }
  
  # display in shaded picture (if 3d)
  if(is.matrix(result) & any(!is.na(result))){
    image.plot(as.numeric(rownames(result)), as.numeric(colnames(result)), result, xlab=names(dimnames(result))[1], ylab=names(dimnames(result))[2], main=outputParameter, col=rev(heat.colors(100)), sub="", graphics.reset=F)
    
    # write simulation numbers into grid
   
      text(as.numeric(rep(rownames(selectedSims), ncol(selectedSims))), as.numeric(rep(colnames(selectedSims), each=nrow(selectedSims))), selectedSims, cex=.8)

    
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