library(fields)

x <- rep(c(1,0,0,0,0,0), 10)
data <- matrix(rep(x, length(x)), ncol=length(x))

#
for(i in 1:nonperiodicBC$parameter$nSteps){
  FFT2d(nonperiodicBC$rasters$vegetation[[i]][-1,-1])
  title(main=paste("timestep", i))
  Sys.sleep(1)
}

punkt <- 40
gesamt <- 400
space <- (gesamt-punkt)/2
data <- matrix(c(rep(rep(0, gesamt),space), rep(c(rep(0,space), rep(1,punkt), rep(0,space)),punkt), rep(rep(0,gesamt),space)), ncol=gesamt)


image.plot(data, asp=1)

FFT2d <- function(data){
  fft_data <- data

  for (i in 1:ncol(data)){
    fft_data[,i] <- fft(data[,i])
  }
  
  #image(Mod(fft_data))
  
  for (i in 1:nrow(data)){
    fft_data[i,] <- fft(fft_data[i,])  
  }
  
  #image.plot(Mod(fft_data), asp=1, nlevel=2000)
  
  # center spectrum:
  cols <- ncol(fft_data)
  rows <- nrow(fft_data)
  fft_centered <- Mod(fft_data)[c((rows/2+1):rows, 1:(rows/2)),c((cols/2+1):cols, 1:(cols/2))]
  image.plot(fft_centered, asp=1, nlevel=2000)
  
  #add 
}


n <- length(fft_data)
sort(Mod(fft_data), partial=n-1)[n-1]

which(Mod(fft_data) == sort(Mod(fft_data), F)[n-0:1])

plot(Mod(fft_data)[1,]
     