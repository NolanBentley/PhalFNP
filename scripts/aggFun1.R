aggFun1 <- function(currFile,byLen=NULL,setParentLen=T){
  library(data.table)
  x<-fread(currFile)
  if(is.null(byLen)){
    currDiffs <- diff(x$V2)
    byLen     <- as.numeric(names(which.max(table(currDiffs[currDiffs!=1]))))
  }
  if(setParentLen){byLen <<- byLen}
  x$interval <- floor((x$V2-x$V2[1])/byLen)+1
  currAg <- aggregate(x$V3,by=list(x$V1,x$interval),function(x){c(mean(x),median(x),sd(x),length(na.omit(x)))})
  outDf <- data.frame(
    win = byLen,
    ind = basename(currFile),
    chr = currAg$Group.1,
    int = currAg$Group.2,
    avg = currAg$x[1],
    med = currAg$x[2],
    sd  = currAg$x[3],
    n   = currAg$x[4]
  )
}
