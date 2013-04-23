
getClusters <- function(grid, threshold){
  # for the function m and n need to be derived inside the function
  clusters <- data.frame(cluster=numeric(), row=numeric(), startcol=numeric(), endcol=numeric(), size=numeric())
  cc <- 0 #cluster count
  n <- ncol(grid)
  m <- nrow(grid)
  
  for(j in 1:m){
    first=T
    for(i in 1:n){
      
      if(grid[j,i]>threshold){ 
        if(first==T){
          cc=cc+1
          clusters[cc,]$cluster <- cc
          clusters[cc,]$startcol <- i
          clusters[cc,]$endcol <- i
          clusters[cc,]$row <- j
          clusters[cc,]$size <- 1
          first=F
        }else{
          clusters[cc,]$endcol <- i
          clusters[cc,]$size <- clusters[cc,]$size + 1
        }
      }else{
        first=T
      }
      
    }
    
  }
  
  for(j in 2:m){
    thisrow <- which(clusters$row==j)
    prevrow <- which(clusters$row==j-1)
    for(a in thisrow){
      for(b in prevrow){
        if(clusters[b,]$startcol <= clusters[a,]$endcol +1  & clusters[b,]$endcol >= clusters[a,]$startcol -1 ){
          clusters[ which(clusters$cluster==clusters[b,]$cluster), ]$cluster <- clusters[a,]$cluster
          
        }
      }
    }
  }
  return(clusters)
}



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