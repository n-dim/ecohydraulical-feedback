source("InputCombination.R")

library(foreach)
#install.packages("doParallel")
library(doParallel)
registerDoParallel(12)


#---- read parameters from file ----
pars <- read.table("/media/Data/eco-hyd/exampleParameters.txt", sep="=", comment.char="!", fill=T, allowEscapes=T, encoding="UTF-8", strip.white=T, as.is=T)

parlist <- t(pars[,2])
colnames(parlist) <- pars[,1]
parlist <- as.list(parlist[1,])

parlist$title <- NULL

save(parlist, file="../example simulation run/exampleParameters.RData")

#---- load parameters ----

load(file="../example simulation run/exampleParameters.RData")

#---- change parameters ----

parlist$m <- parlist$n <- 70
parlist$roughness <- 10^c(-10:13)


#---- write parameter cascade ----

folder <- write.EcoHyd.Input(parlist, threads=12)

#---- run simulation ----

system.time({
foreach(file = list.files(file.path(folder, "input")), .combine=) %dopar% {
  system(paste("../model/ecohydModel.out", file.path(folder, "input",file),  folder) )
}
})


source("CSVtoRData.r")
system.time({
  foreach(file= list.files(folder, pattern="*_inputParameter.txt"), .combine=) %dopar% {
    CSVtoRData(file.path(folder, file))
  }
  
})