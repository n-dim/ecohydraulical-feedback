##---- extractWaveLength ----
extractWaveLength <- function (bild, plot=F, plotSpec=F, smoothing=0) {
  
  suppressPackageStartupMessages({
    library(fields)
    library("spatstat")  
  })
  
    
  Freq <- fft(bild)
  Freq[1] <- 0
  
  # center spectrum:
  cols <- ncol(Freq)
  rows <- nrow(Freq)
  
  rowindex <- floor((-rows/2+1):(rows/2))
  colindex <- floor((-cols/2+1):(cols/2))
  
  fft_centered <- Freq[rowindex%%rows +1, colindex%%cols +1]
  colnames(fft_centered) <- colindex
  rownames(fft_centered) <- rowindex
  
  if(smoothing!=0){
    mod_fft_centered <- as.matrix(blur(as.im(Mod(fft_centered)), smoothing))
  }else{
    mod_fft_centered <- Mod(fft_centered)
  }
  
  # index matrix (using complex numbers) to calculate positions inside the 2d-plane:
  indexMatrix <- matrix(complex(
    real = matrix(rep(colindex, rows), nrow=rows, byrow=F), 
    imaginary = matrix(rep(rowindex, cols), nrow=rows, byrow=T)
  ), nrow=rows)
  
  # show 3d plot:
  if(plotSpec) {library(rgl); persp3d(rowindex/rows,colindex/rows,mod_fft_centered, col="gray")}


  #persp(rowindex/rows,colindex/rows,mod_fft_centered, col="gray", border="gray")
  
  # get the maximum power peak:
  max <- which.max(mod_fft_centered) # maximum power value
  MaxPos <- indexMatrix[max] # position of maximum power value
  rowlength <- Re(MaxPos) # x-position of maximum power position
  collength <- Im(MaxPos) # y-position of maximum power position
  
  phaseshift <- Arg(fft_centered[max]) %% (2*pi)
  #cat("phaseshift:" ,  phaseshift)
  
  orientation <- Arg(MaxPos) %% (2*pi)
  wavelength <- rows/Mod(MaxPos)
  #cat(" wavelength= ", wavelength, "\n")
  
  pointer <- complex(modulus=wavelength, argument=orientation)
  origin <- complex(modulus=rows/2-phaseshift*wavelength/(2*pi), argument=orientation) 

  if(plot){
    image(1:rows, 1:cols, bild, asp=1, col=c("white", "#AAAAAA"), cex.main=1, font.main=1, xlab="m", ylab="n")
    mtext("Image", side=3, line=0)
    mtext(paste0("blue dashed lines (if visible): pattern with wavelength = ", round(wavelength, 3), " pixel and orientation = ", round(orientation, 3), " rad"), side=3, outer=T, padj=3)
    
    #arrows(Re(origin), Im(origin), Re(origin)+ Re(pointer), Im(origin)+Im(pointer))
    
    #---- print the line overlay ----
    lines <- -(abs(collength)+abs(rowlength)):(abs(collength)+abs(rowlength))
    if(length(lines)<60){ # if there are too many lines: don't print
      for(i in lines){
        if(is.infinite(cols/collength)){
          abline(v=i*rows/rowlength-(phaseshift+0+0)*rows/rowlength/(2*pi), col="blue", lty=2)
        }else{
          if(is.infinite(rows/rowlength)){
            abline(h=i*cols/collength-(phaseshift+0+0)*cols/collength/(2*pi), col="blue", lty=2)
          }else{
            abline(a=i*cols/collength-(phaseshift+0+0)*cols/collength/(2*pi), b=-rowlength/collength, col="blue", lty=2)   
          }
        }      
      } 
    }
    image(rowindex, colindex, mod_fft_centered, asp=1, col=gray(100:0/100), xlab="p", ylab="q")
    mtext("Spectrogram", side=3, line=0)
    #if(plot) drape.plot(rowindex/rows,colindex/rows,mod_fft_centered, border=NA, ticktype="detailed", d=2, r=1, phi=30, cex.axis=0.6, cex.lab=0.6, zlab="Amplitude", theta=0, shade=0.1)
  }
  return(list(wavelength=wavelength, orientation=orientation, phaseshift=phaseshift))
}