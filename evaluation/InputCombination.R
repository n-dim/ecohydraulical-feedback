## ---- write.ecoHydInput ----

# example:
file <- "/media/Data/eco-hyd/sim2014-01-06"
data <- list(n=10, m=10, p=4, l=c(3,4,5), t=rep(c(6,7), each=3))

data <- list(a=c(1,4,5), b=c(2,3,20), c=c(10,11,12))

loop <- function(data, Envir, Counter){
  if(Counter$threadSplit==1) Counter$threadSplit = 10000000
  
  for(i in 1:length(data[[1]])){
      assign(names(data[1]), data[[1]][i], envir=Envir)
      #print(paste(names(data[1]),"=", data[[1]][i]))
      if(!length(data)<=1){
        loop(tail(data, -1), Envir, Counter)
      }else{
          objects <- ls(pos=Envir)
          
          if(Counter$simRunNo %% Counter$threadSplit == 1){
           
              Counter$simFile <- paste0(Counter$file, "_No", Counter$simRunNo, "ff.txt")
              write("!eco-hydraulical simulation input file", file=Counter$simFile)

          }
          
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
    Counter$file <- file.path(simFolder,paste0(simName))
    
    combinations <- 1
    for(i in 1:length(data)){
      combinations <- combinations * length(data[[i]])
    }
    Counter$threadSplit <- ceiling(combinations/threads)
       
    loop(data, Envir, Counter)
    
    rm(Envir)
    rm(Counter)
  } 
  
}


#---- read parameters from file ----
pars <- read.table("/media/Data/eco-hyd/exampleParameters.txt", sep="=", comment.char="!", fill=T, allowEscapes=T, encoding="UTF-8", strip.white=T, as.is=T)

parlist <- t(pars[,2])
colnames(parlist) <- pars[,1]
parlist <- as.list(parlist[1,])
as.numeric(parlist,)
nums <- which(!is.na(as.numeric(parlist)))
parlist[nums] <- as.numeric(parlist[nums])

parlist$title <- NULL

#---- change parameters ----

parlist$roughness <- 10^c(-1:1)


#---- write parameter cascade ----

write.EcoHyd.Input(parlist, threads=3)

