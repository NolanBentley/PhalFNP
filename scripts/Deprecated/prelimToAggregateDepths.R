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

#Remove missing values and regions with missing measurements
aggDf <- aggDf[!(is.na(aggDf$avg)|is.na(aggDf$sd)),]
aggDf <- aggDf[aggDf$n>=(0.5*max(aggDf$n)),] #This means this method will not work for homozygous deletions

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
p2 <- ggplot(outDf[abs(outDf$median)>3&grepl("Chr",outDf$chr),],
             aes((max+min)/2000000,median,color=max-min+1))+
  geom_line(color="black")+
  geom_point()+
  scale_color_viridis_c()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  labs(title="All",x="Mbp",y="Median of 5x avg coverage across 5Kbp most 40Kbp")

htmlwidgets::saveWidget(plotly::ggplotly(p1),
                        file=paste0("data_ignored/secondary/depthImages/",aggDf2_rollMed$Group.1[i],".html"))
aggDf$avg_FullNorm_Trimmed[aggDf$regMad>10] <- NA

aggDf$ind_chr <- paste0(aggDf$ind,"_",aggDf$chr)
aggDf_rolled<-aggregate(aggDf$avg_FullNorm_Trimmed,by=list(aggDf$ind_chr),zoo::rollmedian,k=5)

>>>>>>> 167e1b75c5dcdf8a9625f691cb9698defb92affc


#Final region mad
finaRegMad    <- aggregate(aggDf$avg_FullNorm_Trimmed,by=list(aggDf$interval),mad   )
hist(finaRegMad$x,1000);abline(v = c(0.25,2),col="red")
aggDf$finalMad <-finaRegMad$x[match(aggDf$interval,finaRegMad$Group.1)]
aggDf$finalMadLogic <- !(aggDf$finalMad<0.25|aggDf$finalMad>2)
hist(aggDf$finalMad[aggDf$finalMadLogic],1000);abline(v = c(0.25,2),col="red")

library(ggplot2)
aggDf<-aggDf[order(aggDf$ind,aggDf$chr,aggDf$minInt),]
aggDfPlotted <- aggDf[grepl("Chr",aggDf$chr),]
linesToCheck <- c('phal_FIL30_567_M2','phal_FIL30_622_M2','phal_FIL30_846_M2','phal_FIL30_610_M2','phal_FIL30_643_M2','phal_FIL30_83_M2-2')
p1 <- ggplot(aggDfPlotted,
             aes((minInt+maxInt)/2,avg_FullNorm_Trimmed,color=ind,group=ind))+
  geom_line()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  guides(color="none",group="none")+
  scale_y_continuous(limits = c(-15,15))
p2 <- ggplot(aggDfPlotted[aggDfPlotted$finalMadLogic,],
             aes((minInt+maxInt)/2,avg_FullNorm_Trimmed,color=ind,group=ind))+
  geom_line()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  guides(color="none",group="none")+
  scale_y_continuous(limits = c(-15,15))

ggsave(filename = "./data/ploidyPlot5.png",p1,
       width = 20, height = 15, dpi = 600)
print(5)
ggsave(filename = "./data/ploidyPlot6.png",p2,
       width = 20, height = 15, dpi = 600)
print(6)

write.csv(aggDf,file = "./data_ignored/secondary/aggregateDepths4.csv")
 
hist(aggDfPlotted$avg_FullNorm_Trimmed,100000,xlim = c(-1,20),ylim=c(0,1000))

hist(aggDf$avg_FullNorm_Trimmed,100000,xlim=c(-20,40),ylim=c(0,1000))

#Subset to gt 4 normalized
aggDfGt4 <- aggDf[abs(aggDf$avg_FullNorm_Trimmed)>4,]

#Drop regions  with less than 3 consecuttive
aggDfGt4 <- aggDfGt4[order(aggDfGt4$ind,aggDfGt4$chr,aggDfGt4$minInt),]
View(aggDfGt4)

aggDf

p3 <- ggplot(aggDfGt4,
             aes((minInt+maxInt)/2,avg_FullNorm_Trimmed,color=ind,group=ind))+
  geom_line()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  guides(color="none",group="none")+
  scale_y_continuous(limits = c(-15,15))


