#Set variables
opt <- list() #For storing options
opt$wd <- "~/Experiments/PhalFNP/"
opt$prelimFile <- "data_ignored/secondary/prelimAggregateDepths_20240731_0022.csv"
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
aggDf$chr <- gsub("-.*","",aggDf$interval_og)
#aggDf$interval <- paste0(aggDf$interval_og,"-",gsub(" ",0,format(aggDf$minInt)))
aggDf$ind_chr <- paste0(aggDf$ind,"_",aggDf$chr)
aggDf$ind_int <- paste0(aggDf$ind,"_",aggDf$interval)
aggDf <- aggDf[order(aggDf$ind,aggDf$chr,aggDf$minInt),]
aggDf$order <- 1:nrow(aggDf)

#Check for duplicates
aggDf$ind_int_duped <- aggDf$ind_int%in%(aggDf$ind_int[duplicated(aggDf$ind_int)])
aggDf <- aggDf[!duplicated(aggDf$ind_int),]

#Remove errors
aggDf <- aggDf[!(is.na(aggDf$avg)|is.na(aggDf$sd)),]

#Remove intervals with fewer than half samples characterized more than half the length
aggDf$nOverMin <- aggDf$n>=(0.5*max(aggDf$n))
cuL$interval_n <- aggregate(aggDf$interval[aggDf$nOverMin],by=list(aggDf$interval[aggDf$nOverMin]),length)
cuL$nProp      <- cuL$interval_n$x/max(cuL$interval_n$x)
hist(cuL$nProp,1000)
cuL$keptIntervals <- cuL$interval_n$Group.1[cuL$nProp<0.9]
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
missComboDf<-data.frame(file = paste0(gsub("zzzz.*","",names(missingCombos2)),"__",uniInt[missingCombos2],"int"))
missComboDf$interval<-gsub("(.*)(-.of.)","\\1",uniInt[missingCombos2])
missComboDf$ind <- gsub("__.*","",basename(missComboDf$file))
missComboDf$chr <- gsub("-.*","",missComboDf$interval)
missComboDf$interval <- paste0(missComboDf$interval,"-",gsub(" ",0,format(missComboDf$minInt)))
aggDf$ind_chr <- paste0(aggDf$ind,"_",aggDf$chr)
aggDf$ind_int <- paste0(aggDf$ind,"_",aggDf$interval)
aggDf <- aggDf[order(aggDf$ind,aggDf$chr,aggDf$minInt),]
aggDf$order <- 1:nrow(aggDf)

missComboDf<-prepareAggDf(missComboDf)
missComboDf$avg <- 0
missComboDf$med <- 0
missComboDf$sd  <- NA
missComboDf$n   <- 0

## Add in certain values
aggInt          <- aggregate(aggDf$n     ,list(aggDf$interval),median)
aggInt$minInt   <- aggregate(aggDf$minInt,list(aggDf$interval),min   )$x
aggInt$maxInt   <- aggregate(aggDf$maxInt,list(aggDf$interval),max   )$x
missComboDf$rowInAggInt <- match(missComboDf$interval,aggInt$Group.1)
missComboDf$minInt <- aggInt$minInt[missComboDf$rowInAggInt]
missComboDf$maxInt <- aggInt$maxInt[missComboDf$rowInAggInt]

## Add in the fake data
if(!all((missComboDf$ind%in%aggDf$ind)&(missComboDf$chr%in%aggDf$chr))){
  stop("Something went wrong!")
}else{
  aggDf2 <- merge(x = aggDf,y = missComboDf, all= T)
}

#Add in peak values
aggDf$ogSamplePeakMean <- peakFinding(aggDf)

# Remove intervals with low overall n-values
int_n <- aggregate(aggDf2$n,by=list(aggDf2$interval),median)
hist(int_n$x,1000,ylim=c(0,200)); abline(v = 4500,col="red")
aggDf2$IntervalNLogic <- aggDf2$interval%in%int_n$Group.1[int_n$x>4500]
aggDf3 <- aggDf2[aggDf2$IntervalNLogic,]


#### Stopped here I think?

# Roll the means
?zoo::rollmedian
aggDf         <- aggDf[order(aggDf$chr,aggDf$minInt,aggDf$ind),]
aggRollMedian <- aggregate(aggDf$avg,by=list(aggDf$ind_chr),zoo::rollmedian,k=5,fill="extend")

for(i in 1:length(aggRollMedian)){
  if(i==1){outDf <- NULL}
  currLine <- aggRollMedian[[i]]
  for(j in 1:length(currLine$x)){
    outDf<-rbind(outDf,data.frame(i,j,sample=names(avgDepthMat_rolled)[i],chr=currLine$Group.1[j],median=currLine$x[[j]]))
  }
  print(i)
}

outDf_meds <- aggregate(outDf$median,by=list(outDf$sample),median)
outDf$sample_median <- outDf_meds$x[match(outDf$sample,outDf_meds$Group.1)]
outDf$normN <- outDf$median*2/outDf$sample_median
hist(outDf$normN,1000,ylim=c(0,1000))

homoDel <- outDf[outDf$normN<0.5,]

homoDel

########## old 
'
#Remove missing values and regions with missing measurements
#aggDf <- aggDf[!(is.na(aggDf$avg)|is.na(aggDf$sd)),]
#aggDf <- aggDf[aggDf$n>=(0.5*max(aggDf$n)),] #This means this method will not work for homozygous deletions

#Remove regions with strongly skewed distributions
aggDf$meanMinusMedian <- (aggDf$avg - aggDf$med)
aggDiffMed <- aggregate(abs(aggDf$meanMinusMedian),by = list(aggDf$interval),median)
diffNorm <- (aggDiffMed$x-median(aggDiffMed$x))/mad(aggDiffMed$x)
hist(diffNorm,100000,xlim=c(-20,20))
aggDf <- aggDf[aggDf$interval%in%(aggDiffMed$Group.1[diffNorm<20]),]


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
'
