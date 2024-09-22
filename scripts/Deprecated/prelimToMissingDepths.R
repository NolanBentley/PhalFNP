#Set variables
opt <- list() #For storing options
opt$wd <- "~/Experiments/PhalFNP/"
opt$prelimFile <- "data_ignored/secondary/prelimAggregateDepths_20240805_1547.csv"
opt$imageFile  <- "data_ignored/secondary/afterWindow.rimage"
opt$winDfFile  <- "data_ignored/secondary/windowedNs.csv"
opt$samValFile <- "data_ignored/secondary/aggDf_SampleValues.csv"
opt$kVal       <- 7
opt$sourceFile <- "./scripts/functions/peakFinding.R"

#Setup environment
setwd(opt$wd)
library(zoo)
source(opt$sourceFile)
cuL <- list() #For storing analyses

#Load data
aggDf <- read.csv(opt$prelimFile)
aggDf$ind <- gsub("__.*","",basename(aggDf$file))
aggDf$interval_og <- aggDf$interval
aggDf$chr <- gsub("-.*","",aggDf$interval_og)
aggDf$step <- as.numeric(gsub("of.*","",gsub(".*_|int\\.depth\\.out$","",aggDf$file)))
aggDf$interval <- paste0(aggDf$interval_og,"_",gsub(" ","0",format(aggDf$step)))
aggDf$ind_chr <- paste0(aggDf$ind,"_",aggDf$chr)
aggDf$ind_int <- paste0(aggDf$ind,"_",aggDf$interval)
aggDf <- aggDf[order(aggDf$ind,aggDf$chr,aggDf$minInt),]
aggDf$order <- 1:nrow(aggDf)

#Check for duplicates
aggDf$ind_int_duped <- aggDf$ind_int%in%(aggDf$ind_int[duplicated(aggDf$ind_int)])
sum(aggDf$ind_int_duped)

#Remove errors
sum((is.na(aggDf$avg)|is.na(aggDf$sd)))
aggDf <- aggDf[!(is.na(aggDf$avg)|is.na(aggDf$sd)),]

#Remove intervals with fewer than half samples characterized more than half the length
aggDf$nOverMin <- aggDf$n>=(0.5*max(aggDf$n))
cuL$interval_n <- aggregate(aggDf$interval[aggDf$nOverMin],by=list(aggDf$interval[aggDf$nOverMin]),length)
cuL$nProp      <- cuL$interval_n$x/max(cuL$interval_n$x)
hist(cuL$nProp,1000,ylim=c(0,100))
cuL$keptIntervals <- cuL$interval_n$Group.1[cuL$nProp>0.9]
cuL$intervalLogic <- aggDf$interval%in%cuL$keptIntervals; print(mean(cuL$intervalLogic))
aggDf <- aggDf[cuL$intervalLogic,]

#Add in missing intervals
## Generate missing combos
uniInt <- unique(aggDf$interval)
names(uniInt) <- uniInt
missingCombos <- aggregate(aggDf$interval,by=list(aggDf$ind),function(x,y){which(!y%in%x)},y=uniInt)
missingCombosList <- missingCombos$x
names(missingCombosList)<-paste0(missingCombos$Group.1,"zzzz")
missingCombos2 <- unlist(missingCombosList)

## Built data to add
aggDf$status<-"Observed"
simDf <- expand.grid(ind=unique(aggDf$ind),interval=unique(aggDf$interval))
simDf$avg <- 0
simDf$med <- 0
simDf$step <- gsub(".*_","",simDf$interval)
simDf$chr  <- gsub("-.*","",simDf$interval)
simDf$ind_chr <- paste0(simDf$ind,"_",simDf$chr)
simDf$ind_int <- paste0(simDf$ind,"_",simDf$interval)
simDf$order <- 1:nrow(simDf)
simDf$status <- "Sim"
simDf$X <- NA
simDf$sd <- NA
simDf$n <- 0
simDf$minInt <- NA
simDf$maxInt <- NA
simDf$file <- NA
simDf$interval_og <- NA
simDf$ind_int_duped <- NA
simDf$nOverMin <- NA
simDf <- simDf[!simDf$ind_int%in%aggDf$ind_int,]

## Add in the fake data
if(!all(colnames(simDf)%in%colnames(aggDf))){
  stop("Check colnames!")
}else{
  aggDf2 <- merge(x = aggDf,y = simDf, all= T)
}

#####Add analysis values where needed####
#Find sample peaks
aggDf$ogSamplePeakMean <- peakFinding(aggDf)
samDf <- data.frame(ind=unique(aggDf$ind))
samDf$ogSamplePeakMean  <- aggDf$ogSamplePeakMean[match(samDf$ind,aggDf$ind)]
aggDf2$ogSamplePeakMean <- samDf$ogSamplePeakMean[match(aggDf2$ind,samDf$ind)]

#Calculate CN
aggDf2$cn <- aggDf2$avg*2/aggDf2$ogSamplePeakMean

#Find interval min and max pos
intDf     <- aggregate(aggDf$minInt,by=list(aggDf$interval),min)
intDf$max <- aggregate(aggDf$maxInt,by=list(aggDf$interval),max)$x
needMin <- is.na(aggDf2$minInt)
aggDf2$minInt[needMin]<- intDf$x  [match(aggDf2$interval[needMin],intDf$Group.1)]
aggDf2$maxInt[needMin]<- intDf$max[match(aggDf2$interval[needMin],intDf$Group.1)]

#Remove intervals with low overall n-values
int_n   <- aggregate(aggDf2$n  ,by=list(aggDf2$interval),median)
hist(int_n$x ,1000,ylim=c(0,200)); abline(v = max(int_n$x)*0.9,col="red")
aggDf2$nLogic <- aggDf2$interval%in%(int_n$Group.1[int_n$x>=max(int_n$x)*0.9])
aggDf3 <- aggDf2[aggDf2$nLogic,]

#Remove intervals with unusual median CNs
#int_cn  <- aggregate(aggDf2$cn[aggDf2$nLogic],by=list(aggDf2$interval[aggDf2$nLogic]) ,median)
#hist(int_cn$x,10000,xlim=c(0,10)); abline(v = c(1.5,2.5),col="red")
#aggDf2$cnLogic  <- aggDf2$interval%in%(int_n$Group.1[int_cn$x>=1.5&int_cn$x<=2.5])
#aggDf3 <- aggDf2[aggDf2$cnLogic,]

#Report values
c(dims=dim(aggDf3),
  minInt=min(table(aggDf3$interval)),maxInt=max(table(aggDf3$interval)),
  minInd=min(table(aggDf3$ind))     ,maxInd=max(table(aggDf3$ind))
)

#Remove loci with median n_og away from diploid by more than half a copy number 
## Remove ID_Chr combos with different distributions to handle anueploidies
cuL$aggCnMean <- aggregate(abs(aggDf3$cn-2),by=list(aggDf3$ind_chr),mean)
cuL$aggCnMean$chr <- as.numeric(gsub(".*_|Chr","",cuL$aggCnMean$Group.1))
hist(cuL$aggCnMean$x,100000,ylim=c(0,100),xlim=c(0,10))
cuL$perChrIndRank <- aggregate(cuL$aggCnMean$x,by=list(cuL$aggCnMean$chr),rank,ties.method="random")
cuL$aggCnMean$rank<-NA
for(i in 1:nrow(cuL$perChrIndRank)){
  cuL$aggCnMean$rank[cuL$aggCnMean$chr==cuL$perChrIndRank$Group.1[i]] <- cuL$perChrIndRank$x[[i]]
}
cuL$LowOffDiploid_ind_chr<- cuL$aggCnMean$Group.1[cuL$aggCnMean$rank<quantile(cuL$aggCnMean$rank,0.8)]
aggDf3$offDipLogic <- aggDf3$ind_chr%in%cuL$LowOffDiploid_ind_chr

## Calculate and subset
cuL$aggCnIntMedian <- aggregate(aggDf3$cn[aggDf3$offDipLogic],by=list(aggDf3$interval[aggDf3$offDipLogic]),median)
hist(abs(cuL$aggCnIntMedian$x-2),100000,xlim = c(0,10),ylim=c(0,100));abline(v = 0.5,col="red")
cuL$offDiploidLogic <- aggDf3$interval%in%(cuL$aggCnIntMedian$Group.1[abs(cuL$aggCnIntMedian$x-2)<0.5])
mean(cuL$offDiploidLogic)
aggDf3 <- aggDf3[cuL$offDiploidLogic,]

#Recalculate peaks
aggDf3$newSamplePeakMean <- peakFinding(aggDf3)
aggDf3$cn <- aggDf3$avg*2/aggDf3$newSamplePeakMean

#Save the windowed depths
aggDf3 <- aggDf3[order(aggDf3$ind,aggDf3$chr,aggDf3$interval),]
winL <- list()
winL$aggWin_mean  <-aggregate(aggDf3$cn    ,by=list(aggDf3$ind_chr),zoo::rollmean  ,k=opt$kVal)
winL$aggWin_median<-aggregate(aggDf3$cn    ,by=list(aggDf3$ind_chr),zoo::rollmedian,k=opt$kVal)
winL$aggWin_min   <-aggregate(aggDf3$minInt,by=list(aggDf3$ind_chr),function(x,k){zoo::rollapply(x,width=k,FUN=min)},k=opt$kVal)
winL$aggWin_max   <-aggregate(aggDf3$maxInt,by=list(aggDf3$ind_chr),zoo::rollmax   ,k=opt$kVal)

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

#Add discriptors to the window file
hist(diff(aggDf3$min),1000)
winName <- paste0(
  opt$kVal,"x",
  median(aggDf3$max-aggDf3$min+1)/1000,"Kbx",
  median(diff(aggDf3$min))/1000,"Kb_",
  median(winDf$max - winDf$min+1)/1000,"Kb"
)
winDf$winName <- winName
aggDf3$winName <- winName

#Save data
save.image(file = opt$imageFile)
write.csv(winDf,paste0(gsub("\\.csv$","",opt$winDfFile),"_",winName,".csv"))
write.csv(aggDf3[match(unique(aggDf3$ind),aggDf3$ind),],paste0(gsub("\\.csv$","",opt$samValFile),"_",winName,".csv"))

