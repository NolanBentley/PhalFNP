wd <- "~/Experiments/PhalFNP/"
prelimFile <- "data_ignored/secondary/prelimAggregateDepths.csv"

#Load data
setwd(wd)
aggDf <- read.csv(prelimFile)
#aggDf$ind <- gsub("__.*","",basename(aggDf$X))

#Remove missing values and regions with missing measurements
aggDf <- aggDf[!(is.na(aggDf$avg)|is.na(aggDf$sd)),]
aggDf <- aggDf[aggDf$n>=(0.5*max(aggDf$n)),]

#Remove regions with strongly skewed distributions
#aggDf$interval <- aggDf$int
#aggDf$meanMinusMedian <- aggDf$avg - aggDf$med
#aggDiffMed <- aggregate(abs(aggDf$meanMinusMedian),by = list(aggDf$interval),median)
#diffNorm <- (aggDiffMed$x-median(aggDiffMed$x))/mad(aggDiffMed$x)
#hist(diffNorm,10000,xlim=c(-20,20))
#aggDf <- aggDf[aggDf$interval%in%(aggDiffMed$Group.1[diffNorm<20]),]


#Begin transformations
aggMedian <- aggregate(aggDf$avg,by=list(aggDf$ind),median)
aggMad    <- aggregate(aggDf$avg,by=list(aggDf$ind),mad   )

aggDf$sampleMedian<-aggMedian$x[match(aggDf$ind,aggMedian$Group.1)]
aggDf$sampleMad   <-aggMad$x[match(aggDf$ind,aggMad$Group.1)]

aggDf$avg_samNorm <- (aggDf$avg-aggDf$sampleMedian)/aggDf$sampleMad

aggDf$region<- paste0(aggDf$chr,"_",aggDf$interval)
regMedian <- aggregate(aggDf$avg_samNorm,by=list(aggDf$region),median)
regMad    <- aggregate(aggDf$avg_samNorm,by=list(aggDf$region),mad   )
aggDf$regMedian<-regMedian$x[match(aggDf$region,regMedian$Group.1)]
aggDf$regMad   <-regMad$x[match(aggDf$region,regMad$Group.1)]

aggDf$avg_FullNorm <- (aggDf$avg_samNorm-aggDf$regMedian)
aggDf$avg_FullNorm_Trimmed <- aggDf$avg_FullNorm

hist(regMad$x,1000,xlim=c(0,100),ylim=c(0,100));abline(col="red",v = 10)
aggDf$avg_FullNorm_Trimmed[aggDf$regMad>10] <- NA
hist(aggDf$avg_FullNorm_Trimmed,1000,ylim=c(0,100))

#Final region mad
finaRegMad    <- aggregate(aggDf$avg_FullNorm_Trimmed,by=list(aggDf$region),mad   )
hist(finaRegMad$x,1000);abline(v = c(0.25,2),col="red")
aggDf$finalMad <-finaRegMad$x[match(aggDf$region,finaRegMad$Group.1)]
aggDf$finalMadLogic <- !(aggDf$finalMad<0.25|aggDf$finalMad>2)
hist(aggDf$finalMad[aggDf$finalMadLogic],1000);abline(v = c(0.25,2),col="red")

library(ggplot2)
#aggDf$chr <- gsub("-.*","",aggDf$interval)
#aggDf<-aggDf[order(aggDf$ind,aggDf$chr,aggDf$minInt),]
aggDfPlotted <- aggDf[grepl("Chr",aggDf$chr),]
p1 <- ggplot(aggDfPlotted,
             aes(int,avg_FullNorm_Trimmed,color=ind,group=ind))+
  geom_line()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  guides(color="none",group="none")+
  scale_y_continuous(limits = c(-15,15))
p2 <- ggplot(aggDfPlotted[aggDfPlotted$finalMadLogic,],
             aes(int,avg_FullNorm_Trimmed,color=ind,group=ind))+
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

highAgg <- aggDfPlotted[aggDfPlotted$avg_FullNorm_Trimmed>10,]


paste0(gsub("__.*","",head(names(sort(table(highAgg$ind),decreasing = T)))),collapse = "','")


