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

#Save the aggregated / transformed depths
winDf <- data.frame(ind_chr = as.vector(sapply(aggWin_mean$Group.1,rep,times=length(aggWin_mean$x[[1]]))))
winDf$i <- as.vector(sapply(1:nrow(aggWin_mean),rep,times=length(aggWin_mean$x[[1]])))
winDf$id  <- gsub("(^.*)_(.*$)","\\1",winDf$ind_chr)
winDf$chr <- gsub("(^.*)_(.*$)","\\2",winDf$ind_chr)
winDf$mean  <- do.call("c",aggWin_mean$x  )
winDf$median<- do.call("c",aggWin_median$x)
winDf$min   <- do.call("c",aggWin_min$x   )
winDf$max   <- do.call("c",aggWin_max$x   )

####################################
######### Begin Old script #########
####################################

#Begin transformations
## Calculate sample distribution
aggMedian <- aggregate(aggDf$avg,by=list(aggDf$ind),median)
aggMad    <- aggregate(aggDf$avg,by=list(aggDf$ind),mad   )

aggDf$sampleMedian<-aggMedian$x[match(aggDf$ind,aggMedian$Group.1)]
aggDf$sampleMad   <-aggMad$x[match(aggDf$ind,aggMad$Group.1)]

## Normalize by robust sample distribution
aggDf$avg_samNorm <- (aggDf$avg-aggDf$sampleMedian)/aggDf$sampleMad

## Calculate interval distribution
regMedian <- aggregate(aggDf$avg_samNorm,by=list(aggDf$interval),median)
regMad    <- aggregate(aggDf$avg_samNorm,by=list(aggDf$interval),mad   )
aggDf$regMedian<-regMedian$x[match(aggDf$interval,regMedian$Group.1)]
aggDf$regMad   <-regMad$x[match(aggDf$interval,regMad$Group.1)]

## Center by interval mean
aggDf$avg_FullNorm <- (aggDf$avg_samNorm-aggDf$regMedian)
aggDf$avg_FullNorm_Trimmed <- aggDf$avg_FullNorm
aggDf$avg_FullNorm_Trimmed[aggDf$regMad>10] <- NA

#Roll acros windows
aggDf2<-aggDf[!is.na(aggDf$avg_FullNorm_Trimmed),]
aggDf2_rollMea<-aggregate(aggDf2$avg_FullNorm_Trimmed,by=list(aggDf2$ind_chr),zoo::rollmean,k=5)
aggDf2_rollMed<-aggregate(aggDf2$avg_FullNorm_Trimmed,by=list(aggDf2$ind_chr),zoo::rollmedian,k=5)
aggDf2_rollMin<-aggregate(aggDf2$minInt,by=list(aggDf2$ind_chr),FUN=function(x){zoo::rollapply(x,width=5,FUN=min)})
aggDf2_rollMax<-aggregate(aggDf2$maxInt,by=list(aggDf2$ind_chr),zoo::rollmax,k=5)

#install.packages(c("plotly","htmlwidgets","ggplot2"))

dir.create("data_ignored/secondary/depthImages/")
library(ggplot2)
library(plotly)
library(htmlwidgets)

minValue <- min(unlist(lapply(aggDf2_rollMed$x,min,na.rm=T)))
maxValue <- max(unlist(lapply(aggDf2_rollMed$x,max,na.rm=T)))
maxPos   <- max(unlist(lapply(aggDf2_rollMax$x,max,na.rm=T)))
for(i in 1:nrow(aggDf2_rollMea)){
  if(i==1){
    outDf <- NULL
  }
  if(is.null(aggDf2_rollMed$x[[i]][1])){next}
  currDf <- data.frame(
    i = i,
    id = gsub("(^.*)_(.*$)","\\1",aggDf2_rollMed$Group.1[i]),
    chr = gsub("(^.*)_(.*$)","\\2",aggDf2_rollMed$Group.1[i]),
    median = aggDf2_rollMed$x[[i]],
    mean   = aggDf2_rollMea$x[[i]],
    min    = aggDf2_rollMin$x[[i]],
    max    = aggDf2_rollMax$x[[i]]
  )
  outDf <- rbind(outDf,currDf)
  if(i%%1000==1){print(i)}
}

write.csv(outDf,"data_ignored/secondary/windowedDepths5x.csv")
outDf2 <- outDf[abs(outDf$median)>3&grepl("Chr",outDf$chr),]
p2 <- ggplot(outDf2,
             aes((max+min)/2000000,median,color=id))+
  geom_point()+
  guides(color="none")+
  theme_bw()+
  facet_wrap(~chr,ncol = 3)+
  labs(title="All",x="Mbp",y="Median of 5x avg coverage across 5Kbp most 40Kbp")

htmlwidgets::saveWidget(plotly::ggplotly(p2),
                        file=paste0("data_ignored/secondary/above3.html"))

uniChr<-sort(unique(outDf2$chr))
for(i in 1:length(uniChr)){
  p2 <- ggplot(outDf2[outDf2$chr==uniChr[i],],
               aes((max+min)/2000000,median,color=id,fill=max-min+1))+
    geom_point()+
    guides(color="none")+
    theme_bw()+
    facet_wrap(~chr,ncol = 3)+
    labs(title="All",x="Mbp",y="Median of 5x avg coverage across 5Kbp most 40Kbp")
  
  htmlwidgets::saveWidget(plotly::ggplotly(p2),
                          file=paste0("data_ignored/secondary/above3_",uniChr[i],".html"))
  print(i)
}


for(i in unique(outDf2$i)){ 
  p1 <- ggplot(currDf,aes((max+min)/2000000,median,color=max-min+1))+
    geom_line(color="black")+
    geom_point()+
    scale_color_viridis_c()+
    theme_bw()+
    coord_cartesian(xlim=c(0,maxPos),ylim=c(minValue,maxValue))+
    labs(title=aggDf2_rollMed$Group.1[i],x="Mbp",y="Median of 5x avg coverage across 5Kbp most 40Kbp")
  
  htmlwidgets::saveWidget(plotly::ggplotly(p1),
                          file=paste0("data_ignored/secondary/depthImages/",aggDf2_rollMed$Group.1[i],".html"))
  print(i)
}
