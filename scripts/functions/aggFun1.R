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


aggFunTest <- function(currFile,byLen=NULL,setParentLen=T){
  fileSuffix <- gsub(".*_(.of.)int\\.depth\\.out","\\1",test)
  options(datatable.numThreads=1)
  library(data.table)
  x <- fread(currFile)
  if(is.null(byLen)){
    currDiffs <- diff(x$V2)
    byLen     <- as.numeric(names(which.max(table(currDiffs[currDiffs!=1]))))
  }
  if(setParentLen){byLen <<- byLen}
  x$interval <- paste0(fileSuffix,x$V1,"-",gsub(" ",0,format(floor((x$V2-x$V2[1])/byLen)+1)))
  setkey(x,interval)
  
  dtAgg <- x[, .(avg = mean(V3),med = median(V3),sd = sd(V3),n = .N, minInt = min(V2),maxInt = max(V2)), by = interval]
  return(dtAgg)
}




