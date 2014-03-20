#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- source functions ----
source("parameterSpaceFunctions.R")

#-----  -----------

parlist$description <- NULL
parlist$simErosion <- "F"
parlist$roughness <- 0
parlist$dx <- 5
parlist$pa <- 0.2 * 1000
parlist$storEmerge <- 0.048*1000*parlist$dx
parlist$ts <- 1-0.97
parlist$np <- NULL
parlist$n <- NULL
parlist$m <- 32
parlist$kc <- c(0.2,0.4,0.6, 0.8, 1, 2)
parlist$kf <- c(0.2, 0.4, 0.6, 0.8, 1, 2)
parlist$K0 <- 10^seq(-2,1,by=0.2)
parlist$Kmax <- NULL
parlist$Kincrease <- 10^seq(-2,0,by=0.2)

#--- only upper right part of the graphic ----

parlist$kf <- 2
parlist$kc <- 2

################################################
parameterSpace <- calcParameterSpace(parlist)

if(F){
  source("simulate_ParameterSpace.R")
}
