peakFinding <- function(x){
  indOgMedians <- aggregate(x$avg,by=list(x$ind),median)
  indOgMean    <- aggregate(x$avg,by=list(x$ind),mean  )
  indOgMad     <- aggregate(x$avg,by=list(x$ind),mad)
  uniInd <- indOgMedians$Group.1
  if(any(indOgMedians$x<20)){warning("Low og medians present!")}
  x$ogSamplePeakMean <- NA
  plottedIs <- NULL
  plottedIs <- c(plottedIs,highSkew  = which.max(abs(indOgMedians$x-indOgMean$x)))
  plottedIs <- c(plottedIs,lowMedian = which.min(indOgMedians$x))
  plottedIs <- c(plottedIs,highMedian= which.max(indOgMedians$x))
  plottedIs <- c(plottedIs,highMean  = which.max(indOgMean$x))
  plottedIs <- c(plottedIs,lowMad    = which.min(indOgMad$x))
  plottedIs <- c(plottedIs,highMad   = which.max(indOgMad$x))
  try(dev.off(),silent = T);par(mfrow = c(2, 3))
  for(i in 1:length(uniInd)){
    currRows  <- which(x$ind==uniInd[i])
    currAvg   <- x$avg[currRows]
    currAvg   <- currAvg[currAvg<(indOgMedians$x[i]+5*indOgMad$x[i])&currAvg>(indOgMedians$x[i]-5*indOgMad$x[i])]
    currDense <- density(currAvg)
    x$ogSamplePeakMean[currRows] <- currDense$x[which.max(currDense$y)]
    if(i%%10==1){print(i)}
    #Visualize odd samples
    if(i %in% plottedIs){
      plot(currDense,col="red",
           xlim=c(indOgMedians$x[i]-15*indOgMad$x[i],indOgMedians$x[i]+15*indOgMad$x[i]),lwd=3,
           main = paste0(i,": ",names(plottedIs)[plottedIs==i],"_",uniInd[i],collapse = "\n"))
      lines(density(x$avg[currRows]),col="black",lwd=2,lty=2)
      abline(v = c(indOgMedians$x[i],indOgMean$x[i]),col=c("green","purple"),lwd=1,lty=3)
    }
  }
  return(x$ogSamplePeakMean)
}