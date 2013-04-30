coverRatio <- function(data){
  require(fields)
  colsVeg <- two.colors(n=10, start="yellow", end="#003300", middle="green")
  
  with(c(data$parameter, data$rasters), {
  
    coverRatio <- NULL
    mn = (n-1)*(m-1)
    for (i in 1:nSteps){
       coverRatio[i] <- length(which(vegetation[[i]][2:m,2:n]>0))/mn
    }
    
    marold <- par()$mar
    par(mar=c(5,4,4,9))
    
    plot(coverRatio, type="l", xlab="time [y]", ylab="cover ratio [-]", main=paste("simulation run", title), 
         las=1, ylim=c(0,max(coverRatio)))
    Median <- median(coverRatio)
    abline(h=Median, col="red")
    axis(side=2, at=Median, labels=round(Median, digits=2), las=1,col.ticks="red", col="red")
    
    coverRatio1 <- list(NULL)
    Temp <- NULL
    coverRatio1[[1]]  <- rep(0,nSteps)
    coverRatio1$sum  <- rep(0,nSteps)
    
    for(j in 1:9){
      for (i in 1:nSteps){
        Temp[i] <- length(which(vegetation[[i]][2:m,2:n]==j))/mn
      }
      coverRatio1[[j+1]] <- Temp + coverRatio1[[j]]
      plot <- c( coverRatio1[[j+1]] , rev(coverRatio1[[j]]) )
      polygon(c(1:nSteps, nSteps:1), plot, col=colsVeg[j+1], lwd=0.5)
    }
    legend(par()$usr[2]*1.01, mean(par()$usr[3:4]), fill=colsVeg, legend=paste("veg dens", 1:9), xpd=T, yjust=0.5)
    par(mar=marold)
  })
  return(median_cover_ratio=Median)
}