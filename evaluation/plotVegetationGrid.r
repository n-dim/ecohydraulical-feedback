plotVegetationGrid <- function(data, time=1:data$parameter$nSteps){
  require(fields)
  attach(data$parameter)
  
  colsVeg <- c("white", two.colors(n=10, start="yellow", end="green4", middle="green"))
  
  for (i in time){
    image.plot(1:n, 1:m, data$rasters$vegetation[[i]], asp=n/m, col=colsVeg, nlevel=10, 
               main=paste0('simulation run "', title, '", timestep ', i), zlim=c(0,10))
    
  }
  detach(data$parameter)
  
}

#to plot all:
#plotVegetationGrid(Test) 

#to plot only one:
#plotVegetationGrid(Test, 1)

#to plot multiple:
#plotVegetationGrid(Test, c(1,50,100,150,200))