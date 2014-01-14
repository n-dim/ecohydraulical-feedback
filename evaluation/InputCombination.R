## ---- write.ecoHydInput ----

loop <- function(data, Envir, Counter){
  for(i in 1:length(data[[1]])){
      assign(names(data[1]), data[[1]][i], envir=Envir)

      if(!length(data)<=1){
        loop(tail(data, -1), Envir, Counter)
      }else{
          objects <- ls(pos=Envir)

          Counter$simFile <- paste0(Counter$file, "_No", Counter$simRunNo, ".txt")
          write("!eco-hydraulical simulation input file", file=Counter$simFile)
          write(paste("\ntitle =", Counter$simName, " No.", get("simRunNo", Counter), "\n!---------------------\n"), file=Counter$simFile, append=T)

          for(obj in objects){
            write(paste(obj, "=", get(obj,pos=Envir)), file=Counter$simFile, append=T)
          }
          local( simRunNo <- simRunNo + 1, Counter)
      }
  }  
  
}

write.EcoHyd.Input <- function(data, folder= "/media/Data/eco-hyd", simName= "Testrun", threads=1){
  # threads defines the number of parallel computations
  
  if(any(duplicated(names(data)))){
    stop("duplicated entries in input")
  }else{
    
    Envir <- new.env()
    Counter <- new.env()
    Counter$simRunNo <- 1
    Counter$simName <- simName
  
    simFolder <- file.path(folder, paste0("simRun_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S_"), simName))
    dir.create(simFolder)
    dir.create(file.path(simFolder, "input"))
    Counter$file <- file.path(simFolder, "input",paste0(simName))
    
    combinations <- 1
    for(i in 1:length(data)){
      combinations <- combinations * length(data[[i]])
    }
    Counter$threadSplit <- max(round(combinations/threads),1)
       
    loop(data, Envir, Counter)
    
    rm(Envir)
    rm(Counter)
    return(simFolder)
  } 
  
}
