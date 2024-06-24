aggFun1 <- function(currFile,byLen=NULL,setParentLen=T){
  x <- as.data.frame(data.table::fread(currFile))
  if(is.null(byLen)){
    currDiffs <- diff(x$V2)
    byLen     <- as.numeric(names(which.max(table(currDiffs[currDiffs!=1]))))
  }
  if(setParentLen){byLen <<- byLen}
  x$interval <- paste0(x$V1,"-",gsub(" ",0,format(floor((x$V2-x$V2[1])/byLen)+1)))
  currAg  <- aggregate(x$V3,by=list(x$interval,x$V1),function(y){c(mean(y),median(y),sd(y),length(na.omit(y)))})
  currAg2 <- aggregate(x$V2,by=list(x$interval,x$V1),function(y){c(min(y),max(y))})
  outDf <- data.frame(
    win = byLen,
    ind = basename(currFile),
    chr = currAg$Group.2,
    int = currAg$Group.1,
    avg = currAg$x[1],
    med = currAg$x[2],
    sd  = currAg$x[3],
    n   = currAg$x[4],
    intMin = currAg2$x[1],
    intMax = currAg2$x[2]
  )
  outDf<-outDf[order(outDf$chr,outDf$intMin,outDf$int,outDf$ind),]
  return(outDf)
}

aggFun2 <- function(currFile,byLen=NULL,setParentLen=T){
  options(datatable.numThreads=1)
  library(data.table)
  x <- fread(currFile)
  if(is.null(byLen)){
    currDiffs <- diff(x$V2)
    byLen     <- as.numeric(names(which.max(table(currDiffs[currDiffs!=1]))))
  }
  if(setParentLen){byLen <<- byLen}
  x$interval <- paste0(x$V1,"-",gsub(" ",0,format(floor((x$V2-x$V2[1])/byLen)+1)))
  setkey(x,interval)
  
  dtAgg <- x[, .(avg = mean(V3),med = median(V3),sd = sd(V3),n = .N, minInt = min(V2),maxInt = max(V2)), by = interval]
  return(dtAgg)
}



