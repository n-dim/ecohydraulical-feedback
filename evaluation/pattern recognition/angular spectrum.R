## ---- angular_spectrum ----
get_angular_spectrum <- function (indexMatrix, powerMatrix, breaks, radiusBreaks , plotDiagram=T, plotImage=T, thresholdFrac=0.01) {
  # indexMatrix ... complex matrix containing the index numbers for the power matrix
  # powerMatrix ... matrix containing the power (modulus of the fft)
  # breaks ... vector of radial breaks to bin the result
  # radiusBreaks ... if only the angular spectrum of a certain interval of radii shall be calculated
  # plotDiagram ... plot a diagram
  # plotImage ... plot an image of the radial segmentation and their weights
  # thresholdFrac ... Fraction of the interval 1/length(breaks)...1 to set the threshold line
  

  library("fields")
  
  angularMatrix <- Arg(indexMatrix) %% (2*pi)
  breaks  <- breaks %% (2*pi)
  #powerMatrix <- powerMatrix - mean(powerMatrix) + 1
  
  angular_diagram <- NULL
  angularPlot  <- angularMatrix*NA
  mids <- head(breaks,-1) + (tail(breaks, -1) - head(breaks, -1))/2 # mids between the breaks
  
  for(i in 2:length(breaks)){
    if(breaks[i] > breaks[i-1]){
      inAngle <- which(
        angularMatrix  <= breaks[i] & 
        angularMatrix > breaks[i-1] &
        Mod(indexMatrix) <= max(radiusBreaks) &
        Mod(indexMatrix) > min(radiusBreaks)
      )    
    }else{
      inAngle <- which(
        angularMatrix  <= breaks[i]  |
        angularMatrix > breaks[i-1] &
        Mod(indexMatrix) <= max(radiusBreaks) &
        Mod(indexMatrix) > min(radiusBreaks)
      ) 
    }
    
    angularPlot[inAngle] <- angular_diagram[i-1] <- max(powerMatrix[inAngle])
  }
  angular_diagram <- angular_diagram/sum(angular_diagram)
  
  whiteNoise <- 1/length(breaks)
  threshold <- whiteNoise + (1-whiteNoise)*thresholdFrac
  overthreshold <- which(angular_diagram > threshold)
    
  #--- plots ---
  if(plotDiagram) {
    plot(mids, angular_diagram, type="l", main="Angular Spectrum", xlab="angle [rad]", ylab="power", las=1, ylim=range(angular_diagram)*c(1,1.4))
    
    abline(h=whiteNoise, lty=3)
    abline(h=threshold, lty=2)
    points(mids[overthreshold], angular_diagram[overthreshold])
    
    legend("top", lty=c(3,2), c("hypothetic white noise distribution", paste0(thresholdFrac*100, "% deviation from white noise distribution")), box.lty=0, bg="white", cex=.8)
    #threshold <- quantile(angular_diagram, probs=c(.8,.85,.9), na.rm=T)
    #abline(h=threshold, lty=2)
    #axis(side=4, at=threshold, labels=names(threshold), las=1)
  #points(breaks[which.max(angular_diagram)], max(angular_diagram, na.rm=T), col="red")
  }
  if(plotImage) image.plot(Im(indexMatrix)[1,], Re(indexMatrix)[,1], angularPlot, main="Angular Spectrum", asp=1, col=gray(90:0/100), xlab="x", ylab="y", legend.lab="power")
  
  
  
  return(list(breaks=breaks, mids=mids, power=angular_diagram, thresholdFrac=thresholdFrac, threshold=threshold, overthreshold=overthreshold))
}