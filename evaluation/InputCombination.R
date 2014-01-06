## ---- write.ecoHydInput ----

# example:
file <- "/media/Data/eco-hyd/sim2014-01-06"
data <- list(n=10, m=10, p=4, l=c(3,4,5), t=rep(c(6,7), each=3))

data <- list(a=c(1,4,5), b=c(2,3,20), c=c(10,11,12))

loop <- function(data, simName){
  for(i in 1:length(data[[1]])){
      assign(names(data[1]), data[[1]][i], envir=Envir)
      #print(paste(names(data[1]),"=", data[[1]][i]))
      if(!length(data)<=1){
        loop(tail(data, -1), simName)
      }else{
          objects <- ls(pos=Envir)
        
          write(paste("\n!----", simName, "- No.", get("simRunNo", Counter), "\n"), file=file, append=T)
          for(obj in objects){
            write(paste(obj, "=", get(obj,pos=Envir)), file=file, append=T)
          }
          local( simRunNo <- simRunNo + 1, Counter)
      }
  }  
  
}

write.EcoHyd.Input <- function(data, folder= "/media/Data/eco-hyd", simName= "Testrun"){
  if(any(duplicated(names(data)))){
    stop("duplicated entries in input")
  }else{
    
    Envir <- new.env()
    Counter <- new.env()
    assign("simRunNo", 1, envir=Counter)
  
    simFolder <- file.path(folder,paste0("simRun_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S_"), simName))
    dir.create(simFolder)
    file <- file.path(simFolder,paste0(simName, ".txt"))
    
    
    write("!--- input variation ---", file=file)
    
    loop(data, simName)
    
    rm(Envir)
    rm(Counter)
  } 
  
}

write.EcoHyd.Input(data)


n=10
m=10
p=4
for(l in c(3,4,5)){
  
}


write.ecoHydInput <- function(data, file){
  
  for(i in 1:nrow(data)){
    write.table(t(data[i,]), file, sep=" = ", quote=F, append=(!i==1), col.names=F)
    write("\n!--------------------------------\n", file=file, append=T)
  }
  
}

write.ecoHydInput(data, file)

#read.table(file=file, sep="=", header=F)
