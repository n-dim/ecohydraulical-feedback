coverRatioComparison <- function(path, name){
  source("coverRatio.r")
  
  files <- dir(path, pattern="grids.RData")
  
  dataEnv <-  new.env()
  for(i in files){
    load(file.path(path,i), dataEnv)
  }
  
  pdf(file=file.path(path,paste("coverRatioComparison", name, ".pdf")), width=11.96, height=8.27)
    datasets <- ls(envir=dataEnv)
  
    for ( i in datasets){
      coverRatio(get(i, envir=dataEnv))
    } 
  dev.off()
  
}