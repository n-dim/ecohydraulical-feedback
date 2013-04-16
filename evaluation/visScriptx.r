library(fields)

readCSV <- function(file=NA,nrow=64,ncol=64,tsteps=c(200)) {
  # Choose File with input parameters to read in first
  if(is.na(file)) file <- file.choose()
  nskip <- tsteps+(tsteps-1)*nrow
  data <- array(NA,dim=c(nrow,ncol,length(tsteps)))
  for ( i in 1:length(tsteps)){
        dat <- scan(file,sep=";",strip.white=T,skip=nskip[i],nlines=nrow)
        data[,,i] <- dat[!is.na(dat)]
  }
  data
  }
  

tsteps <- c(200)
nrows <- 64
ncols <- 64
  veg1 <- readCSV(file="/home/hydro/Desktop/output/01-vegetation.csv")
  veg2 <- readCSV(file="/home/hydro/Desktop/output/02-vegetation.csv")
  veg3 <- readCSV(file="/home/hydro/Desktop/output/03-vegetation.csv")
  veg4 <- readCSV(file="/home/hydro/Desktop/output/04-vegetation.csv")
  
  
pdf(file="/home/hydro/Desktop/output/exampleEco.pdf")
colsVeg <- two.colors(n=9, start="yellow", end="green4", middle="green")
par(mfrow=c(2,2))

    image(veg1[,,1],asp=1,zlim=c(1,9),col=colsVeg)
    image(veg2[,,1],asp=1,zlim=c(1,9),col=colsVeg)
    image(veg3[,,1],asp=1,zlim=c(1,9),col=colsVeg)
    image(veg4[,,1],asp=1,zlim=c(1,9),col=colsVeg)
dev.off()

