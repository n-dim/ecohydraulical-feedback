##---- entropy ----
get_entropy <- function(x){
  N <- length(x)
  Hmax <- -log2(1/N)
  H <- -sum(x*log2(x), na.rm=T)
  return(H/Hmax)
}
  
