cntCall<-function(x){sum(!is.na(x))}
rollmin<-function(x,k,fill){zoo::rollapply(data=x,width=k,fill=fill,FUN=min    )}
rollN  <-function(x,k,fill){zoo::rollapply(data=x,width=k,fill=fill,FUN=cntCall)}

badLocusFilter <- function(indDepth,badLoci){
  for(i in 1:nrow(badLoci)){
    currLocus <- badLoci[i,]
    lociToRemove <-
      indDepth$chrN == currLocus$chrN&
      indDepth$start<=(currLocus$mid+badLocusSpacing)&
      indDepth$end  >=(currLocus$mid-badLocusSpacing)
    indDepth <- indDepth[!lociToRemove,]
  }
  return(indDepth)
}

prepareDepths <- function(i,currFile,currId,k,lociToRemove,chrToNumPattern = "Chr|scaffold_",filterToChr = T,chrPattern ="Chr",filterOnBadLoci = F){
  #Read and clean
  indDepth   <- read.delim(currFile) 
  if(!(nrow(indDepth)>0&ncol(indDepth)==6)){stop("Input not correct dimensions")}
  colnames(indDepth)<- c("chr","start","end","mappable","counts","CN")
  
  #Prep data
  if(filterToChr){indDepth <- indDepth[grep(chrPattern,indDepth$chr),]}
  indDepth$chrN <- as.numeric(gsub(chrToNumPattern,"",indDepth$chr))
  indDepth$i    <- i
  indDepth$file <- currFile
  indDepth$id   <- currId
  indDepth$order<- 1:nrow(indDepth)
  indDepth$mid  <- (indDepth$start+indDepth$end)/2
  
  #Filter the bad loci if appropriate
  if(filterOnBadLoci){
    indDepth <- badLocusFilter(indDepth,lociToRemove)
  }
  
  #Make a rolled analysis
  indDepth <- cbind(indDepth,data.frame(
    winMeanChr=zoo::rollmean(  indDepth$chrN  ,k = k,fill = T),
    winMinPos =     rollmin(   indDepth$start ,k = k,fill = T),
    winMaxPos =zoo::rollmax(   indDepth$end   ,k = k,fill = T),
    winMeanCN =zoo::rollmean(  indDepth$CN    ,k = k,fill = T),
    winMedCN  =zoo::rollmedian(indDepth$CN    ,k = k,fill = T),
    winIntN   =       rollN(   indDepth$start ,k = k,fill = T)
  ))
  indDepth$winLogic  <- ((indDepth$winMeanChr%%1)==0)&(indDepth$winIntN==k)
  indDepth$winLen    <-  indDepth$winMaxPos-indDepth$winMinPos+1
  indDepth$winMid    <- (indDepth$winMaxPos+indDepth$winMinPos)/2
  
  #Return data
  return(indDepth)
}

formatId <- function(id,width){
  idPre <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\1",id)
  idNum <- as.numeric(gsub("(.*FIL..)_(\\d+)_(.*)","\\2",id))
  idSuf <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\3",id)
  idNumBuffered <- gsub(" ","0",format(c(10^width,idNum))[-1])
  idForm    <- paste0(idPre,"_",idNumBuffered,"_",idSuf)
  idNum_bin <- gsub(" ","0",paste0("Lines",format(c(10^width,floor(idNum/100)*100))[-1],"-",format(c(10^width,(floor(idNum/100)+1)*100-1))[-1]))
  return(data.frame(id_form=idForm,idNum_bin=idNum_bin))
}

prepareCNVForAggPlot <-function(x,winName){
  x$id_orig  <- x$id
  x$id       <- x$id_form
  x$chr_orig <- x$chr
  x$chr      <- gsub("scaffold.*","scaffold",x$chr_orig)
  x$min      <- x$winMinPos
  x$max      <- x$winMaxPos
  x$IntervalMidpoint_Mbp <- x$mid/1000000
  x$median   <- x$winMedCN
  x$hasDivergentMedian <- x$divLogic
  x$singleLineMultiChr <- paste0(x$id,"_",x$chr,"_",winName,".png")
  x$multiLineSingleChr
  return(x)
}
