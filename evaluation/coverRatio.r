coverRatio <- function(data){
  require(fields)
  colsVeg <- c("black", two.colors(n=10, start="yellow", end="darkgreen", middle="green"))
  
  with(c(data$parameter, data$rasters), {
  
    coverRatio <- NULL
    mn = n*m
    for (i in 1:nSteps){
       coverRatio[i] <- length(which(vegetation[[i]]>0))/mn
    }
    
    marold <- par()$mar
    par(mar=c(5,4,4,9))
    
    plot(coverRatio, type="l", xlab="time [y]", ylab="cover ratio [-]", main=paste("simulation run", title), 
         las=1, ylim=c(0,max(coverRatio)))
    Median <- median(coverRatio)
    abline(h=Median, col="red")
    axis(side=2, at=Median, labels=round(Median, digits=2), las=1,col.ticks="red", col="red")
    
    coverRatio1 <- NULL
    for(j in 1:9){
      for (i in 1:nSteps){
        coverRatio1[i] <- length(which(vegetation[[i]]==j))/mn
      }
      lines(coverRatio1, col=colsVeg[j+1])
    }
    legend(par()$usr[2]*1.01, mean(par()$usr[3:4]), lty=1, col=colsVeg, legend=c('sum', paste("veg dens", 1:9)), xpd=T, yjust=0.5)
    par(mar=marold)
  })
}