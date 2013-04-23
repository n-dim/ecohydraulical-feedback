setwd("/home/nanu/Dokumente/simRuns/2013-04-16/output/")

files <- dir("/home/nanu/Dokumente/simRuns/2013-04-16/output/", pattern=".RData")
for(file in files[1]){
  load(file)
  
}
load('02.RData')
load('01.RData')
Run <- `02`
load('03.RData')
Run <- `03`
load('04.RData')
Run <- `04`

attach(Run$parameter)
attach(Run$rasters)

coverRatio <- NULL
mn = n*m
for (i in 1:nSteps){
   coverRatio[i] <- length(which(vegetation[[i]]>0))/mn
}
plot(coverRatio, type="l", xlab="time [y]", ylab="cover ratio [-]", main=paste("simulation run", title), las=1)
Median <- median(coverRatio)
abline(h=Median, col="red")
axis(side=2, at=Median, labels=round(Median, digits=2), las=1,col.ticks="red", col="red")

detach(Run$rasters)
detach(Run$parameter)

rm(list = ls())
