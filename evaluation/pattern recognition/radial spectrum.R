## ---- radial_spectrum ----
get_radial_spectrum <- function (indexMatrix, powerMatrix, breaks, plotDiagram=T, plotImage=T, thresholdFrac=0.01) {
  # indexMatrix ... complex matrix containing the index numbers for the power matrix
  # powerMatrix ... matrix containing the power (modulus of the fft)
  # breaks ... vector of radial breaks to bin the result
  # plotDiagram ... plot a diagram
  # plotImage ... plot an image of the radial segmentation and their weights
  # thresholdFrac ... Fraction of the interval 1/length(breaks)...1 to set the threshold line
  
  radiusMatrix <- Mod(indexMatrix)
      
  radial_diagram <- NULL
  radiusPlot  <- radiusMatrix*NA
  mids <- head(breaks, -1) + (tail(breaks,-1) - head(breaks, -1))/2
  
  for(i in 2:length(breaks)){
    inRadius <- which(
      radiusMatrix < breaks[i] & 
      radiusMatrix >= breaks[i-1]
    )
    if(length(inRadius)!=0){
      radiusPlot[inRadius] <- radial_diagram[i-1] <- max(powerMatrix[inRadius], na.rm=T)
    }else{
      radiusPlot[inRadius] <- radial_diagram[i-1] <- NA
    }
  }
  radial_diagram <- radial_diagram/sum(radial_diagram, na.rm=T)
  
  whiteNoise <- 1/length(breaks)
  threshold <- whiteNoise + (1-whiteNoise)*thresholdFrac
  overthreshold <- which(radial_diagram > threshold)
  
  #--- plots ---
  if(plotDiagram){
    plot(mids, radial_diagram, type="l", main="Radial Spectrum", xlab="frequency", ylab="power", ylim=range(radial_diagram)*c(1,1.3))
    
    abline(h=whiteNoise, lty=3)
    abline(h=threshold, lty=2)
    points(mids[overthreshold], radial_diagram[overthreshold])
    
    legend("top", lty=c(3,2), c("hypothetic white noise distribution", paste0(thresholdFrac*100, "% deviation from white noise distribution")), box.lty=0, bg="white", cex=.8)
  }
  if(plotImage) image.plot(Im(indexMatrix)[1,], Re(indexMatrix)[,1], radiusPlot, main="Radial Spectrum", asp=1, col=gray(100:0/100), xlab="x", ylab="y", legend.lab="power")
  
  
  return(list(breaks=breaks, mids=mids, power=radial_diagram, thresholdFrac=thresholdFrac, threshold=threshold, overthreshold=overthreshold))
}