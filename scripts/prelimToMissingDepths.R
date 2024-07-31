wd <- "~/Experiments/PhalFNP/"
prelimFile <- "data_ignored/secondary/prelimAggregateDepths_20240731_0022.csv"

#Load data
setwd(wd)
aggDf <- read.csv(prelimFile)
prepareAggDf <- function(x){
  x$ind           <- gsub("__.*","",basename(x$file))
  if(is.null(x$interval_orig)){x$interval_orig <- x$interval}
  x$chr           <- gsub("-.*","",x$interval_orig)
  x$interval      <- paste0(x$interval_orig,"-",gsub(".*([[:digit:]]+of[[:digit:]]+)int.*","\\1",(x$file)))
  x$ind_chr       <- paste0(x$ind,"_",x$chr)
  x$ind_int       <- paste0(x$ind,"_",x$interval)
  return(x)
}
aggDf <- prepareAggDf(aggDf)
aggDf <- aggDf[order(aggDf$chr,aggDf$minInt,aggDf$ind),]
aggDf$order <- 1:nrow(aggDf)

hist(table(aggDf$interval),1000)

#Add in missing intervals
## Generate missing combos
uniInt <- unique(aggDf$interval)
names(uniInt) <- uniInt
missingCombos <- aggregate(aggDf$interval,by=list(aggDf$ind),function(x,y){which(!y%in%x)},y=uniInt)
missingCombosList <- missingCombos$x
names(missingCombosList)<-paste0(missingCombos$Group.1,"zzzz")
missingCombos2 <- unlist(missingCombosList)
## Built data to add
missComboDf<-data.frame(X = paste0(gsub("zzzz.*","",names(missingCombos2)),"__",uniInt[missingCombos2],"int"))
missComboDf$interval<-gsub("(.*)(-.of.)","\\1",uniInt[missingCombos2])
missComboDf<-prepareAggDf(missComboDf)
missComboDf$interval
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
## Checking inputs
if(!all(missComboDf$ind%in%aggDf$ind&missComboDf$chr%in%aggDf$chr)){stop("Something went wrong!")}
## Add in the fake data
aggDf2 <- merge(x = aggDf,y = missComboDf, all= T)


# Remove intervals with low overall n-values
int_n <- aggregate(aggDf2$n,by=list(aggDf2$interval),median)
hist(int_n$x,1000,ylim=c(0,200)); abline(v = 4500,col="red")
aggDf2$IntervalNLogic <- aggDf2$interval%in%int_n$Group.1[int_n$x>4500]
aggDf <- aggDf2[aggDf2$IntervalNLogic,]

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
