#Set variables
opt <- list() #For storing options
opt$wd <- "~/Experiments/PhalFNP/"
opt$prelimFile <- "data_ignored/secondary/prelimAggregateDepths_20240625_1442.csv"
opt$imageFile  <- "data_ignored/secondary/afterWindow.rimage"
opt$winDfFile  <- "data_ignored/secondary/windowedNs5x.csv"
source("~/Experiments/PhalFNP/scripts/peakFinding.R")
cuL <- list() #For storing analyses

#Setup environment
setwd(opt$wd)

#Load data
aggDf <- read.csv(opt$prelimFile)
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
cuL$interval_n <- aggregate(aggDf$interval,by=list(aggDf$interval),length)
cuL$keptIntervals <- cuL$interval_n$Group.1[(cuL$interval_n$x/max(cuL$interval_n$x))>0.5]
cuL$intervalLogic <- aggDf$interval%in%cuL$keptIntervals; print(mean(cuL$intervalLogic))
aggDf <- aggDf[cuL$intervalLogic,]

#Remove regions with strongly skewed distributions
aggDf$pctDiffSkew <- (aggDf$avg - aggDf$med)/(sum(aggDf$avg + aggDf$med)/2)
cuL$aggDiffMed <- aggregate(aggDf$pctDiffSkew,by = list(aggDf$interval),median)
cuL$diffNorm   <- (cuL$aggDiffMed$x-median(cuL$aggDiffMed$x))/mad(cuL$aggDiffMed$x)
hist((cuL$diffNorm),100000,xlim=c(-40,40))
abline(v = c(-20,20),col="red")
cul$intervalLogic2 <- aggDf$interval%in%(cuL$aggDiffMed$Group.1[abs(cuL$diffNorm)<20])
aggDf <- aggDf[cul$intervalLogic2,]

#Find peak depth across each sample
aggDf$ogSamplePeakMean <- peakFinding(aggDf)

#Create ploidy transformed avg 
aggDf$n_og <- aggDf$avg*2/aggDf$ogSamplePeakMean
try(dev.off(),silent = T);hist(aggDf$n_og,100000,xlim=c(0,10),ylim=c(0,1000))

#Remove loci with median n_og away diploid by more than half a copy number 
## Remove IDs combos with different distributions
cuL$aggNMean   <- aggregate(aggDf$n_og,by=list(aggDf$ind_chr),mean)
cuL$aggNMean$chr <- gsub(".*_","",cuL$aggNMean$Group.1)
cuL$aggNMean$absDiffDip <- abs(cuL$aggNMean$x-2)
hist(cuL$aggNMean$absDiffDip,10000)
cuL$aggNOffDipRank <- aggregate(cuL$aggNMean$absDiffDip,by=list(cuL$aggNMean$chr),rank,ties.method="random")
cuL$aggNMean$rank<-NA
for(i in 1:nrow(cuL$aggNOffDipRank)){
  cuL$aggNMean$rank[cuL$aggNMean$chr==cuL$aggNOffDipRank$Group.1[i]] <- cuL$aggNOffDipRank$x[[i]]
}
cuL$LowOffDiploid_ind_chr<- cuL$aggNMean$Group.1[cuL$aggNMean$rank<quantile(cuL$aggNMean$rank,0.8)]
aggDf$aggDfODLogic <- aggDf$ind_chr%in%cuL$LowOffDiploid_ind_chr

## Calculate and subset
cuL$aggNIntMedian <- aggregate(aggDf$n_og[aggDf$aggDfODLogic],by=list(aggDf$interval[aggDf$aggDfODLogic]),median)
cuL$aggNIntMad    <- aggregate(aggDf$n_og[aggDf$aggDfODLogic],by=list(aggDf$interval[aggDf$aggDfODLogic]),mad   )
hist(abs(cuL$aggNIntMedian$x-2),10000)
plot(cuL$aggNIntMedian$x-2,cuL$aggNIntMad$x)
abline(v = 0.5,col="red")
cuL$offDiploidLogic <- aggDf$interval%in%(cuL$aggNIntMedian$Group.1[abs(cuL$aggNIntMedian$x-2)<0.5]); print(mean(cuL$offDiploidLogic))
aggDf <- aggDf[cuL$offDiploidLogic,]

#Recalculate peaks
aggDf$newSamplePeakMean <- peakFinding(aggDf)
aggDf$n <- aggDf$avg*2/aggDf$newSamplePeakMean

#Save the windowed depths
winL <- list()
winL$aggWin_mean  <-aggregate(aggDf$n_og  ,by=list(aggDf$ind_chr),zoo::rollmean  ,k=5)
winL$aggWin_median<-aggregate(aggDf$n_og  ,by=list(aggDf$ind_chr),zoo::rollmedian,k=5)
winL$aggWin_min   <-aggregate(aggDf$minInt,by=list(aggDf$ind_chr),function(x,k){zoo::rollapply(x,width=k,FUN=min)},k=5)
winL$aggWin_max   <-aggregate(aggDf$maxInt,by=list(aggDf$ind_chr),zoo::rollmax   ,k=5)

#Assemble aggregated values into a data.frame
winDf <- data.frame(
  mean  = do.call("c",winL$aggWin_mean$x  ),
  median= do.call("c",winL$aggWin_median$x),
  min   = do.call("c",winL$aggWin_min$x   ),
  max   = do.call("c",winL$aggWin_max$x   )
)
winDf$ind_chr <- unlist(mapply(function(x,y){rep(x,length(y))},x=winL$aggWin_mean$Group.1,y=winL$aggWin_mean$x))
winDf$id  <- gsub("(^.*)_(.*$)","\\1",gsub("_scaffold","",winDf$ind_chr))
winDf$chr <- gsub("(^.*)_(.*$)","\\2",winDf$ind_chr)

#Save data
save.image(file = opt$imageFile)

#Write to a file
write.csv(winDf,opt$winDfFile)
