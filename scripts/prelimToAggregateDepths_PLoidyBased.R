#Setup environment
wd <- "~/Experiments/PhalFNP/"
prelimFile <- "data_ignored/secondary/prelimAggregateDepths_20240625_1442.csv"

#Load data
setwd(wd)
aggDf <- read.csv(prelimFile)
aggDf$ind <- gsub("__.*","",basename(aggDf$X))
aggDf$chr <- gsub("-.*","",aggDf$interval)
aggDf$interval <- paste0(aggDf$interval,"-",gsub(" ",0,format(aggDf$minInt)))
aggDf$ind_chr <- paste0(aggDf$ind,"_",aggDf$chr)
aggDf <- aggDf[order(aggDf$ind,aggDf$chr,aggDf$minInt),]
aggDf$order <- 1:nrow(aggDf)

#Remove observations with fewer than half the bases with reads
aggDf <- aggDf[!(is.na(aggDf$avg)|is.na(aggDf$sd)),]
aggDf <- aggDf[aggDf$n>=(0.5*max(aggDf$n)),] #This means this method will not work for homozygous deletions

#Remove intervals with fewer than half samples characterized
interval_n <- aggregate(aggDf$interval,by=list(aggDf$interval),length)
keptIntervals <- interval_n$Group.1[(interval_n$x/max(interval_n$x))>0.5]
aggDf <- aggDf[aggDf$interval%in%keptIntervals,]

#Remove regions with strongly skewed distributions
aggDf$pctDiffSkew <- (aggDf$avg - aggDf$med)/(sum(aggDf$avg + aggDf$med)/2)
aggDiffMed <- aggregate(aggDf$pctDiffSkew,by = list(aggDf$interval),median)
diffNorm <- (aggDiffMed$x-median(aggDiffMed$x))/mad(aggDiffMed$x)
hist((diffNorm),10000,xlim=c(-40,40))
abline(v = c(-20,20),col="red")
aggDf <- aggDf[aggDf$interval%in%(aggDiffMed$Group.1[abs(diffNorm)<20]),]

#Find peak depth across each sample
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

#Create ploidy transformed avg 
aggDf$ogSamplePeakMean <- peakFinding(aggDf)
aggDf$n_og <- aggDf$avg*2/aggDf$ogSamplePeakMean
try(dev.off(),silent = T);hist(aggDf$n_og,100000,xlim=c(0,10),ylim=c(0,1000))

#Remove loci with median n_og away diploid by more than half a copy number 
## Remove IDs combos with different distributions
aggNMean   <- aggregate(aggDf$n_og,by=list(aggDf$ind_chr),mean)
aggNMean$chr <- gsub(".*_","",aggNMean$Group.1)
aggNMean$absDiffDip <- abs(aggNMean$x-2)
hist(aggNMean$absDiffDip,10000)
aggNOffDipRank <- aggregate(aggNMean$absDiffDip,by=list(aggNMean$chr),rank,ties.method="random")
aggNMean$rank<-NA
for(i in 1:nrow(aggNOffDipRank)){
  aggNMean$rank[aggNMean$chr==aggNOffDipRank$Group.1[i]] <- aggNOffDipRank$x[[i]]
}
LowOffDiploid_ind_chr<- aggNMean$Group.1[aggNMean$rank<quantile(aggNMean$rank,0.8)]
aggDfODLogic <- aggDf$ind_chr%in%LowOffDiploid_ind_chr

## Calculate and subset
aggNIntMedian <- aggregate(aggDf$n_og[aggDfODLogic],by=list(aggDf$interval[aggDfODLogic]),median)
aggNIntMad    <- aggregate(aggDf$n_og[aggDfODLogic],by=list(aggDf$interval[aggDfODLogic]),mad   )
hist(abs(aggNIntMedian$x-2),10000)
plot(aggNIntMedian$x-2,aggNIntMad$x)
abline(v = 0.5,col="red")
aggDf <- aggDf[aggDf$interval%in%(aggNIntMedian$Group.1[abs(aggNIntMedian$x-2)<0.5]),]

#Recalculate peaks
aggDf$newSamplePeakMean <- peakFinding(aggDf)
aggDf$n <- aggDf$avg*2/aggDf$newSamplePeakMean

#Do a window based analysis of consecutive intervals
aggDf$rollMedian<-aggregate(aggDf2$avg_FullNorm_Trimmed,by=list(aggDf2$ind_chr),zoo::rollmean,k=5)

#Save the windowed depths
zoo_rollmin  <- function(x,k){zoo::rollapply(x,width=k,FUN=min)}
aggWin_mean  <-aggregate(aggDf$n_og  ,by=list(aggDf$ind_chr),zoo::rollmean  ,k=5)
aggWin_median<-aggregate(aggDf$n_og  ,by=list(aggDf$ind_chr),zoo::rollmedian,k=5)
aggWin_min   <-aggregate(aggDf$minInt,by=list(aggDf$ind_chr),zoo_rollmin    ,k=5)
aggWin_max   <-aggregate(aggDf$maxInt,by=list(aggDf$ind_chr),zoo::rollmax   ,k=5)

rm(ls()[!ls()%in%c("aggDf","zoo_rollmin","aggWin_mean","aggWin_median","aggWin_min","aggWin_max")])

#Assemble aggregated values into a data.frame
winDf <- data.frame(ind_chr = as.vector(sapply(aggWin_mean$Group.1,rep,times=length(aggWin_mean$x[[1]]))))
winDf$i <- as.vector(sapply(1:nrow(aggWin_mean),rep,times=length(aggWin_mean$x[[1]])))
winDf$id  <- gsub("(^.*)_(.*$)","\\1",winDf$ind_chr)
winDf$chr <- gsub("(^.*)_(.*$)","\\2",winDf$ind_chr)
winDf$mean  <- do.call("c",aggWin_mean$x  )
winDf$median<- do.call("c",aggWin_median$x)
winDf$min   <- do.call("c",aggWin_min$x   )
winDf$max   <- do.call("c",aggWin_max$x   )

#Write to a file
write.csv(winDf,"data_ignored/secondary/windowedNs5x.csv")