## ---- centFFT ----
centFFT <- function(image){
  # calculate Fourier transform
  FFT <- fft(image)
  
  # prepare centering
  cols <- ncol(FFT) 
  rows <- nrow(FFT)
  rowindex <- 1:rows - ceiling(rows/2) 
  colindex <- 1:cols - ceiling(cols/2)
  
  # center FFT
  FFT <- FFT[rowindex %% rows + 1, colindex %% cols + 1]
  rownames(FFT) <- rowindex
  colnames(FFT) <- colindex
  p <- rep(rowindex, cols)
  q <- rep(colindex, each=rows)
  
  # set mean = 0
  FFT["0","0"] <- 0
  
  # get amplitudes
  FFTamplitudes <- Mod(FFT)
  
  # transform to fraction of total variance
  FFTamplitudes <- FFTamplitudes / sum(FFTamplitudes) * 2 
  
  
  # prepare calculation of polar coordinates
  complexIndex <- complex(real=p, imaginary=q)
  
  # calculate radius from complex plane
  radius <- Mod(complexIndex)
  
  # calculate angle from complex plane
  angle <- Arg(complexIndex) %% (2*pi)
  
  return(list(p=rowindex, q=colindex, amplitudes=FFTamplitudes, phaseShift=Arg(FFT), complex=FFT , radius=matrix(radius, rows), angle=matrix(angle, rows), complexIndex=matrix(complexIndex, rows), pMat=matrix(p, rows), qMat=matrix(q,rows)))
}