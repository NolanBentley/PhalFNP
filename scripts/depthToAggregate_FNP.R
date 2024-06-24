#Setup variables
wd       <- "~/Experiments/PhalFNP/"
depthDir <- "./data_ignored/secondary/depths"
outDir   <- "./data/"
prelimFile <- paste0("prelimAggregateDepths_",format(Sys.time(), "%Y%m%d_%H%M"),".csv")
totalSteps <- 5
nCores     <- min(c(10,parallel::detectCores()/2))

#Setup the environment
setwd(wd)
library(parallel)

## Load functions
# Calculate the mean depth per window (byLen)
source('./scripts/aggFun1.R')

#Get data
depthFileList <- lapply(1:totalSteps,function(x){list.files(depthDir,pattern = paste0(x,"of",totalSteps,"int.depth\\.out$"),full.names = T)})
cl <- makeCluster(nCores)
j<-1
for(j in 1:length(depthFileList)){
  if(j==1){aggDf <- NULL}
  byLen <- NULL #Triggers aggFun1 to recalculate each step
  cat("\n\n\n",j,"\n")
  currDepthFiles <- depthFileList[[j]]
  #currDepthFiles<-currDepthFiles[1:40]
  parOut   <- parallel::parSapply(cl,currDepthFiles,aggFun2,simplify = F)
  parOutDf <- do.call("rbind",parOut)
  #View(parOutDf)
  aggDf <- rbind(aggDf,parOut)
}
write.csv(aggDf,file = file.path(outDir,prelimFile))
stopCluster(cl)

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