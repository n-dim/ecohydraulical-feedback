

clusterSizeProgression(data)

grid <- Test$rasters$vegetation[[200]]>0

image(grid)

clusters <- getClusters(grid)

clusterSizes <- getClusterSizes(clusters) 


summary(clusterSizes$size)

with(Test$parameter,
     
     boxplot(clusterSizesStat, main=paste("simulation run", title))
     
)
plotVegetationGrid(Test, time=200)


##############################################

library(fields)
colsVeg <- two.colors(n=9, start="yellow", end="green4", middle="green")

for(i in 1:4){
  
  grid <- Run03$rasters$vegetation[[i]]
  # n <- ncol(grid)
  # m <- nrow(grid)
  # image(1:n, 1:m, grid,zlim=c(1,9),col=colsVeg)
  
  threshold <- 0
  th_grid <- matrix(as.numeric(grid<=threshold), ncol=ncol(grid))
  
  system.time(clusters <- getClusters(grid, threshold))
  
  #pdf(file=Run03$parameter$title, paper="a4")
  image(1:m,1:n,th_grid, col=c("green","white"), main=paste("simulation run ", Run03$parameter$title, "timestep ", i))
  text(labels=clusters$cluster, x=clusters$row, y=clusters$startcol, cex=0.5)
  #dev.off()
  
  agClusters <- aggregate(size ~ cluster, data=clusters, FUN="sum")
  cat('number of clusters: ', nrow(agClusters))
  cat('Cluster sizes: ')
  print(summary(agClusters$size))
  
}