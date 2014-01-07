source("InputCombination.R")

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

parlist$roughness <- 10^c(-1:1)


#---- write parameter cascade ----

folder <- write.EcoHyd.Input(parlist, threads=3)

#---- start simulation ----

for(file in list.files(folder)){
  system(paste("../model/ecohydModel.out", file.path(folder,file),  folder))
}
