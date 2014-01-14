Dimnames <- list(n=c(10,20,30),m=c(10,20,30),A=c("0.1","0.2","0.3"), roughness=c(0.001,0.01,0.1))
Dims <- mapply(length, Dimnames)


parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)

parameterSpace

parameterSpace[1,1,,]

dimnames(parameterSpace)

colnames(parameterSpace[15])
dims[15]
dims(arrayInd(15, c(3,3,3,3), dims))


#parameterSpace[1,,,, drop=T]

for(simNo in 1:length(parameterSpace)){
  pos <- which(parameterSpace==simNo, arr.ind=T)
  parameterSet <- mapply(`[`, dimnames(parameterSpace), pos)
  cat("\ntitle =",simNo, fill=T)
  for(i in 1:length(parameterSet)){
    cat(paste(names(parameterSet[i]), "=", parameterSet[i]), fill=T)
  }
}