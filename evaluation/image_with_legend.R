image.with.legend <- function (..., prettyScale=T, breaks=NA, colorScale="gray.colors") {
  #col <- list(...)$col

  layout(matrix(2:1, nrow=1), widths = c(1,lcm(3.5)))
  
  # draw legend
  if(all(is.na(breaks))){
    if(prettyScale){
      scale <- pretty(c(min(x), max(x)))
    }else{
      scale <- seq(min(x), max(x), length.out=200)
    }
  }else{
    if(length(breaks)==1){
      if(prettyScale){
        scale <- pretty(c(min(x), max(x)), n=breaks)
      }else{
        scale <- seq(min(x), max(x), length.out=breaks+1)
      }
    }else{
       scale <- breaks
    }
  }
    
  tryCatch(expr, error = function(e) e)
  
  col <- rev(unlist(lapply((length(scale)-1), colorScale)))
  
  
  marold <- par()$mar
  par(mar=c(marold[1]+2,0,marold[3]+2,4))
  #image(matrix(scale, nrow=1), axes=F, frame.plot=T)
  plot(c(0,1), range(scale), type="n", axes=F, frame.plot=T, xaxs="i", yaxs="i", ann=F)
  rect(rep(0, length(scale)-1), head(scale, -1) , rep(1, length(scale)-1), tail(scale, -1), col=col, border=NA)
  if(any(is.na(breaks))){
    axis(4, pretty(scale))
  }else{
    axis(4, (scale))
  }
  
  
  box()

  par(mar=marold)
  # draw image
  image(..., breaks=scale, col=col)  
  box()
  
}


# just for testing purpose:
if(FALSE){
  x <- matrix(runif(5^2, 0,3), nrow=5)
  
  image.with.legend(1:5, 1:5, x, prettyScale=T, breaks=NA, colorScale="gray.colors", axes=F)
  
  axis(1, 1:5, letters[1:5])  
}

