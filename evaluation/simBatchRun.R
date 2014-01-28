simBatchRun <- function (parameterSpace, folder, cores=12) {
  #folder ... folder to create the simulation run folder into
  #cores  ... number of threads to use to do multiple simulation runs at once (should match the machine's cores)
  
  # for parallel computing:
  library("foreach")
  library("doParallel")
  registerDoParallel(cores)
  source("coverRatio.r")
  source("readCSV.r")
  
  
  #---- create simulation folder ----
  simFolder <- file.path(folder, paste0("simRun_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))
  dir.create(simFolder)
  
  cat("starting time                       ", format(Sys.time(), "%Y-%m-%d %T"), file=file.path(simFolder,"proceeding.txt"), fill=T)
  
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
    
    # run simulation and write into logfile
    system(paste("../model/ecohydModel.out", file,  simNoFolder, ">", paste0(simNoFolder, "/", simNo, "_log.txt")))
    
    # read into R format
    Data <-  readCSV(file.path(simNoFolder, paste0(simNo, "_inputParameter.txt")))
    
    Data$postp$coverRatio <- coverRatio(Data)
    
    Data$postp$coverRatioMedian <- median(Data$postp$coverRatio)
    
    # save Data
    outFileName <- paste0("sim", simNo)
    assign(outFileName, Data)
    save(list=outFileName, file=file.path(simNoFolder,paste(simNo, "_grids.RData", sep="")))  
    
    cat("finished simulation no", format(simNo, width=3), "of", format(length(parameterSpace),width=3), "at", format(Sys.time(), "%Y-%m-%d %T"), file=file.path(simFolder,"proceeding.txt"), fill=T, append=T)
    
    Data
  }

  cat("end time                            ", format(Sys.time(), "%Y-%m-%d %T"), file=file.path(simFolder,"proceeding.txt"), fill=T, append=T)
  return(simSpace)
}