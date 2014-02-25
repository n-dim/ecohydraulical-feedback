##---- analyseSpectrum ----
#depends on source("../SpektroAnalyse/entropy.R")
analyseSpectrum <- function (image, plotImage=T, plotSpectrum=T, plotRadial=T, plotAngular=T, thresholdFrac=0.03) {
  
  if(plotImage){
    m=1:nrow(image)
    n=1:ncol(image)
    image(m, n, image, col=gray(100:0/100), asp=1, main="Image")
  } 
  FFT <- centFFT(image)
  
  #calculate entropy:
  entropy_2D <- get_entropy(FFT$amplitudes[1:(nrow(FFT$amplitudes)/2),])
  
  if(plotSpectrum){
    with(FFT,{
      image(p,q,amplitudes, col=gray(100:0/100), asp=1, main="Spectrogram")
      title(sub=paste("H/Hmax =", round(entropy_2D,5)), xpd=NA)    
    })
  }
  
  radial_spectrum <- get_radial_spectrum(FFT$complexIndex, FFT$amplitudes, 1:max(FFT$radius), plotDiagram=plotRadial , F, thresholdFrac=thresholdFrac)
  entropy_radial <- get_entropy(radial_spectrum$power)
  if(plotRadial) title(sub=paste("H/Hmax =", round(entropy_radial,3)),xpd=NA)
  
  angleBreaks <- seq(5,185, by=10)/360*2*pi 
  angular_spectrum <- get_angular_spectrum(FFT$complexIndex, FFT$amplitudes, angleBreaks, radiusBreaks=c(0,200), plotDiagram=plotAngular, plotImage=F, thresholdFrac=thresholdFrac)  
  entropy_angular <- get_entropy(angular_spectrum$power)
  if(plotAngular) title(sub=paste("H/Hmax =", round(entropy_angular,3)), xpd=NA)
  
  
  return(list(radial=radial_spectrum, angular=angular_spectrum, entropy_2D=entropy_2D, entropy_radial=entropy_radial, entropy_angular=entropy_angular))
}