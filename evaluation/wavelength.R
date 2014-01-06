#extract wave speed and wave length from image time series

require(fields)

waveSpeed <- function(x,y){ 
  #convert to time series 
  xt <- as.ts(x) 
  yt <- as.ts(y) 
  #calculates cross spectrum 
  xy.spc <- spec.pgram(ts.union(xt, yt),plot=FALSE)   
  
  posPeak <- xy.spc$spec[,1]==max(xy.spc$spec[,1]) 
  wavelength <- 1/xy.spc$freq[posPeak] 
  phs <- xy.spc$phase[posPeak] 
  vel <- phs/2/pi*wavelength 
  return(c(vel,wavelength)) 
} 


col.veg <- two.colors(n=9,start="yellow",middle="Green",end="darkgreen")
nm <- 32
nts <- 400
dx <- 5
dt <- 1

allveg <-  array(dat,dim=c(nm,nm,nf))


# Read File
  
  ############load data
  strm <-file(paste(paste('TStepResults',k,sep=""),'_1.out',sep=""),open = "r")
  readLines(con=strm,n=1)
  veg        <-    array(data=NA,dim=c(nm,nm,nts),dimnames=c("x","y","t"))
  for (i in 1:nts) {
    readLines(con=strm,n=1);
    dat <- scan(file=strm,what=numeric(0),nlines=nm,quiet=TRUE);
    dat <- matrix(data=dat,nrow=nm,ncol=nm*8,byrow=TRUE)
    dat<-array(dat,dim=c(nm,nm,8))
    veg[,,i]        <- dat[,,1]
  }
  close(strm)  
  allveg[,,k] <- veg[,,nts]
  
# new input read:

  load("/home/nanu/simRuns/periodic Boundaries/periodicBC_grids.RData")
  data <- periodicBC
  attach(data$rasters)
  attach(data$parameter)
  veg        <-    array(data=NA,dim=c(m,n,nSteps),dimnames=c("x","y","t"))
  for(i in 1:nSteps){
    veg[,,i] <- vegetation[[i]]
  }
  
############calc wave speeds and wave lengths
  wavelength <- c()
  wavespeed <- c()
  spd <- c()
  wlnth <- c()
  pos <- cbind(1:m,1:n)
  for (i in 2:nSteps) {
    x <- c()
    y <- c()
    for (j in 1:n){
      x <- c(x,veg[pos[j,1],pos[j,2],i-1])
      y <- c(y,veg[pos[j,1],pos[j,2],i])
    }
    res <- waveSpeed(x,y)
    spd <- c(spd,res[1])
    wlnth <- c(wlnth,res[2])
  }
  spd <- spd*(2^0.5)*5
  wlnth <- wlnth*(2^0.5)*5
  wavespeed <- cbind(wavespeed,spd)
  wavelength <- cbind(wavelength,wlnth)
  


###extract means of each

sp.mn <- apply(wavespeed[floor(2*nSteps/3):(nSteps-1),],MARGIN=2,mean,na.rm=TRUE)  # for mean wave length uses only 2/3 of the simulation timeline
sp.mn <- mean(wavespeed[floor(2*nSteps/3):(nSteps-1),])
wl.mn <- apply(wavelength[floor(2*nSteps/3):(nSteps-1),],MARGIN=2,mean,na.rm=TRUE)
wl.mn <- mean(wavelength[floor(2*nSteps/3):(nSteps-1),])

sp.mn.mat <- matrix(sp.mn,nrow=8,byrow=TRUE)
wl.mn.mat <- matrix(wl.mn,nrow=8,byrow=TRUE)

kf <- seq(0.1,1.5,by=0.2)
kc <- kf

contour(x = kc, y = kf, wl.mn.mat,
        nlevels = 5,
        labcex = 1.0, drawlabels = TRUE, method = "simple",
        axes = TRUE,
        col = par("fg"), lty = par("lty"), lwd = par("lwd")
)

contour(x = kc, y = kf, sp.mn.mat,
        nlevels = 5,
        labcex = 1, drawlabels = TRUE, method = "simple",
        axes = TRUE,
        col = par("fg"), lty = par("lty"), lwd = par("lwd")
)


########plot each pattern   
# par(oma=c(0.1,0.1,0.1,0.1))
# par(mar=c(0.05,0.05,0.05,0.05))
# par(mfrow=c(3,4))
# par(ask=TRUE)
# for (i in 1:nf) {
#   image(allveg[,,i],asp=1,xlim=c(0,1),ylim=c(0,1),zlim=c(1,9),axes=FALSE,col=col.veg)
#  }

#non-periodic  manually derived from visual inspection
is.banded <- c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,
               TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,
               TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,
               TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
               TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
               TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
               TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,
               FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)

#periodic manually derived from visual inspection
#is.banded <- c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,
#               TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,
#               TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
#               TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
#               TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
#               TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
#               TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,
#               FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)

is.banded.mat <- matrix(is.banded,nrow=8,byrow=TRUE)  

sp.mn.mat2 <- sp.mn.mat
wl.mn.mat2 <- wl.mn.mat
for ( i in 1:8){
  for (j in 1:8){
    if (is.banded.mat[i,j]==FALSE){
      sp.mn.mat2[i,j]<- NA
      wl.mn.mat2[i,j]<- NA  
    }   
  }}

col.cont <- two.colors(n=3,start="lightgrey",middle="lightgrey",end="white")
postscript(file="/path/waveLengthNonPeriodic2.eps",width=5,height=5)
image(kc,kf, is.banded.mat,,col=col.cont,xlab=expression(tilde(k)[c]),ylab=expression(tilde(k)[f]),cex.lab=1.2)
contour(x = kc, y = kf, wl.mn.mat2,
        levels = c(30,40,50,70,90,110),
        labcex = 1.0, drawlabels = TRUE, method = "simple",
        axes = TTRUE,
        col = par("fg"), lty = par("lty"), lwd = par("lwd")
        ,add=TRUE,xlim=c(0.1,1.5),ylim=c(0.1,1.5))
axis(side=1,lwd=2)
axis(side=2,lwd=2)
dev.off()


postscript(file="/path/waveSpeedNonPeriodic2.eps",width=5,height=5)
image(kc,kf, is.banded.mat,col=col.cont,xlab=expression(tilde(k)[c]),ylab=expression(tilde(k)[f]),cex.lab=1.2)
contour(x = kc, y = kf, sp.mn.mat2,
        levels = c(0,0.5,1.0,2.0,4.0),
        labcex = 1, drawlabels = TRUE, method = "simple",
        axes = TRUE,
        col = par("fg"), lty = par("lty"), lwd = par("lwd")
        ,add=TRUE,xlim=c(0.1,1.5),ylim=c(0.1,1.5))
axis(side=1,lwd=2)
axis(side=2,lwd=2)
dev.off()   


# new wavelength metod:

load("~/simRuns/periodic Boundaries/nonperiodicBC_grids.RData")
attach(nonperiodicBC$rasters)
attach(nonperiodicBC$parameter)

library(doMC) # for multicore calculations
registerDoMC(detectCores()) 

t=nSteps
getWaveLength <- function(grid, dx=1){
  x_spec <- foreach(x=1:nrow(grid), .combine=c) %dopar% {
    per <- spec.pgram(grid[x,], plot=F)
    dominant_frequency <- 1/per$freq[which.max(per$spec)]
  }
  
  wavelength_x <- median(x_spec)
  
  y_spec <- foreach(y=1:ncol(grid), .combine=c) %dopar% {
    per <- spec.pgram(grid[,y], plot=F)
    dominant_frequency <- 1/per$freq[which.max(per$spec)]
  }
  
  wavelength_y <- median(y_spec)
  
  wavelength <- wavelength_x * wavelength_y /sqrt(wavelength_x^2 + wavelength_y^2)
  wavelength * dx
}

wavelengthProgression <- function(data, time=NA){
  with(c(data$rasters, data$parameter),{
    if(is.na(time)) time=1:nSteps
    dx * foreach(t=time, .combine=c) %dopar% getWaveLength(vegetation[[t]])
  })   
}

wavelengthProgression(nonperiodicBC)
plot(wavelengthProgression(nonperiodicBC), type="l", ylab="wavelength [m]", xlab="time [y]")


# Experiments:

x <- rep(c(1,0,0,0,0,0), 10)
data <- matrix(rep(x, length(x)), ncol=length(x))
image(data)
getWaveLength(data)

data <- matrix(rep(x, length(x)-1), ncol=length(x))
image(data)
getWaveLength(data)

data <- matrix(rep(x, length(x)-2), ncol=length(x))
image(data)
getWaveLength(data)

#.-------------------
periodogram(vegetation[[200]][,60]) -> per
1/per$freq[which.max(per$spec)]
