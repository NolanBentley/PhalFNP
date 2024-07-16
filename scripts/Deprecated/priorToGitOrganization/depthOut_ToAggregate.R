#Visualize depths
library(data.table)
depthVec2 <- list.files("~/Experiments/FastNeutron/analysis/depths/",pattern = "depth\\.out$",full.names = T)
#Calculate the mean depth per window (byLen)
aggFun1 <- function(currFile,byLen){
  library(data.table)
  x<-fread(currFile)
  x$interval <- floor(x$V2/byLen)
  currAg <- aggregate(x$V3,by=list(x$V1,x$interval),function(x){c(mean(x),sd(x),length(na.omit(x)))})
  outDf <- data.frame(
    win = byLen,
    ind = basename(currFile),
    chr = currAg$Group.1,
    int = currAg$Group.2,
    avg = currAg$x[,1],
    sd  = currAg$x[,2],
    n   = currAg$x[,3]
  )
}
i<-1
for(i in 1:length(depthVec2)){
  if(i==1){aggDf <- NULL}
  aggDf  <- rbind(aggDf,aggFun1(depthVec2[i],byLen))
  print(i)
}

write.csv(aggDf,file = "~/Experiments/FastNeutron/analysis/prelimAggregateDepths.csv")
#aggDf <- data.table::fread("~/Experiments/FastNeutron/analysis/prelimAggregateDepths.csv")

#Remove missing values abnd regions with missing measurements
aggDf <- aggDf[!(is.na(aggDf$avg)|is.na(aggDf$sd)),]
aggDf <- aggDf[aggDf$n>=quantile(aggDf$n,0.05),]

aggMedian <- aggregate(aggDf$avg,by=list(aggDf$ind),median)
aggMad    <- aggregate(aggDf$avg,by=list(aggDf$ind),mad   )

aggDf$sampleMedian<-aggMedian$x[match(aggDf$ind,aggMedian$Group.1)]
aggDf$sampleMad   <-aggMad$x[match(aggDf$ind,aggMad$Group.1)]

aggDf$avg_samNorm <- (aggDf$avg-aggDf$sampleMedian)/aggDf$sampleMad

aggDf$region<- paste0(aggDf$chr,"_",aggDf$int)
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
p1 <- ggplot(aggDf,
             aes(int,avg_FullNorm_Trimmed,color=ind,group=ind))+
  geom_line()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  guides(color="none",group="none")+
  scale_y_continuous(limits = c(-15,15))
p2 <- ggplot(aggDf[aggDf$finalMadLogic,],
             aes(int,avg_FullNorm_Trimmed,color=ind,group=ind))+
  geom_line()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  guides(color="none",group="none")+
  scale_y_continuous(limits = c(-15,15))

ggsave(filename = "~/Experiments/FastNeutron/analysis/ploidyPlot5.png",p1,
       width = 20, height = 15, dpi = 600)
ggsave(filename = "~/Experiments/FastNeutron/analysis/ploidyPlot6.png",p2,
       width = 20, height = 15, dpi = 600)

write.csv(aggDf,file = "~/Experiments/FastNeutron/analysis/aggregateDepths3.csv")