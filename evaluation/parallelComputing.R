foreach()

library(multicore)
library(microbenchmark)

cl <- makeCluster(4)

system.time({
    for(i in 1:4){
      parallel(microbenchmark( sqrt(1:1000), times=1000000))
    }
    
  collect()
})


system.time({
  for(i in 1:2){
    microbenchmark( sqrt(1:1000), times=1000000)
  }
})

stopCluster(cl)

#    for (i in 1:nSteps){
#       clusters[[i]] <- parallel(getClusters(vegetation[[i]][2:m,2:n]>0), "t")
#       if (i%%cores==0){
#         for(j in (i-cores+1):i){
#           clusters[[j]] <- collect(clusters[[j]])
#         }
#       }
#    }
