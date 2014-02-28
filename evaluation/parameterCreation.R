#---- load parameters ----
load(file="../example simulation run/exampleParameters.RData")

#---- source functions ----
source("parameterSpaceFunctions.R")

#---- change parameters ----

parlist$description <- NULL
parlist$m <- parlist$n <- 60
#parlist$pa <- seq(200, 1100, by=50)
parlist$pa <- c(400,600, 800,1200)
#parlist$pa <- 600
parlist$roughness <- 0.0000001
parlist$Emax <- 876
#parlist$Emax <- seq(200, 2000, by=100)
#parlist$roughness <- 0
#parlist$nSteps <- 100
parlist$BCs <- "-3\t-3\t0\t0"
parlist$simErosion <- "F"
#parlist$np <- 150
parlist$np <- 0 # np is then calculated as pa/4
parlist$nSteps <- 100
#parlist$K0 <- 0.2546103 + seq(0.0, 10, by=1)
#parlist$K0  <- c(0.005, 0.05, 0.5, 1, 2)
parlist$K0 <- seq(0.5, 1, by=0.2)
#parlist$Kmax <- 2.0461027
#parlist$Kmax <- seq(0.5, 6, by=0.5)
parlist$Kmax <- NULL # to calc Kmax as K0 + Kincrease or Kmax * KincFrac
#parlist$Kincrease <- seq(0.1, 2, by=0.2)
parlist$Kincrease <- NULL
#parlist$Kincrease <- c(1, 2, 3, 4)
#parlist$KincFrac <- 1 + 2^c(-2,-1,0,1,2)
parlist$KincFrac <- 1 + 2^c(2,3,4)
parlist$kf <- 1.0
parlist$kc <- 0.6
parlist$useRandomSeed = "T"
#parlist$m <- parlist$n <- 100
#parlist$pa <- seq(200, 1200, by=20)
#parlist$Kmax <- 2.0461027 + seq(0.1, 2, by=0.1)

#--- run no 5----
parlist$KincFrac <- seq(5,15, 2)
parlist$K0 <- 0.5

#--- run no 6 ---
parlist$pa <- c(400,600, 800)
parlist$KincFrac <- c(4,9,15,16,17) 

#--- run no 7 ---
parlist$pa <- 600
parlist$KincFrac <- c(2,3,4,9,15,16,17,18) 
parlist$K0 <- c(0.025, 0.5)

#--- run no 8 ---
parlist$pa <- 600
parlist$KincFrac <- c(2,3,4,9,15,16,17,18) 
parlist$K0 <- seq(0.06, 0.5, by=0.04)

#--- run no 9pre ---
parlist$KincFrac <- c(4,4.5,5,6,9) 
parlist$K0 <- 0.26

#--- run no 9 ---
parlist$KincFrac <- c(4,5,6,9,15,16,17,18,19,20) 
parlist$K0 <- seq(0.02, 0.5, by=0.04)

#--- run with higher m,n for better wavelength detection ----

parlist$m  <- parlist$n  <- 120

#--- run with higher nSteps ----

parlist$nSteps <- 200

#--- check repeatability ----

parlist$K0 <- 0.14
parlist$KincFrac <- 9
parlist$useRandomSeed <- rep(F, 48)