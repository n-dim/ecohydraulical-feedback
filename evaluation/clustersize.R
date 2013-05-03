
getClusters <- function(grid){
  # for the function m and n need to be derived inside the function
  clusters <- data.frame(cluster=numeric(), row=numeric(), startcol=numeric(), endcol=numeric(), size=numeric())
  cc <- 0 #cluster count
  n <- ncol(grid)
  m <- nrow(grid)
  
  for(j in 1:m){
    first=T
    for(i in 1:n){
      
      if(grid[j,i]==1){ 
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

getClusterSizes <- function(clusters){
  return(aggregate(size ~ cluster, data=clusters, FUN="sum"))
} 

clusterSizeProgression <- function(data){
  library(doMC) # multicore calculations
  registerDoMC(detectCores()) 
  
  with(c(data$rasters, data$parameter), {
   
   clusters <- foreach(i=1:nSteps) %dopar% getClusters(vegetation[[i]][2:m,2:n]>0)

   clusterSizes <- list(NULL)
   clusterSizesStat <- as.data.frame(t(rep(0,6)))
   names(clusterSizesStat) <- names(summary(0))
   
   for( i in 1:nSteps){
     clusterSizes[[i]] <- getClusterSizes(clusters[[i]]) 
     clusterSizesStat[i,] <- summary(clusterSizes[[i]]$size) 
   }
   plot(clusterSizesStat$Mean, type="l", ylab="mean cluster size [cells]", xlab="time [y]", main=paste("simulation run", title))
  })
}
