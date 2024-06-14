#Setup variables
wd    <-"~/Experiments/FastNeutron/analysis"
bamPath <- "~/Experiments/FastNeutron/downloading/"
byLen <- 200000
startPos <- 50
regionWidth <- 5000

initDf <- data.frame(
  chr=c('Chr01','Chr02','Chr03','Chr04','Chr05','Chr06','Chr07','Chr08','Chr09'),
  max=c('59908650','57868950','64542600','52912950','59521150','45727350','48048900','44948900','73933400')
)

#Setup packages
packageDir <- file.path(wd,"R")
packageVec <- c("ggplot2","chron","data.table")


#Use paths
setwd(wd)
dir.create("bed"     ,showWarnings = F)
dir.create("depths"  ,showWarnings = F)
dir.create(packageDir,showWarnings = F)
.libPaths(packageDir)
if(!all(packageVec%in%installed.packages())){install.packages(packageVec,lib=packageDir, repos='http://cran.us.r-project.org')}
library("chron")


#Find bam files
bamVec<-list.files(bamPath,pattern = "bam$",full.names = T)

#Remove files being downloaded still
bamVec_infor <- file.info(bamVec)
bamVec <- bamVec[as.chron(bamVec_infor$mtime)< (as.chron(Sys.time())-1/24)]
max(as.chron(file.info(bamVec)$mtime))

#Create filenames
byLenText <- paste0(gsub("\\+","",format(byLen)),"x",regionWidth)
bedName   <- file.path(wd,"bed",paste0("regions_",byLenText,".bed"))
depthVec  <- file.path(wd,"depths",paste0(gsub("\\.bam$","",basename(bamVec)),"_",byLenText,"int.depth.out"))

#Remove files with the files already extant
bamVec <- bamVec[!file.exists(depthVec)]
bamVec_infor <- bamVec_infor[!file.exists(depthVec),]
depthVec <- depthVec[!file.exists(depthVec)]


#Generate bed file
for(i in 1:nrow(initDf)){
  currDf <- data.frame(
    chr=initDf$chr[i],
    min=seq(startPos,initDf[i,2],by=byLen),
    max=seq(startPos,initDf[i,2],by=byLen)+regionWidth
  )
  if(i==1){
    outDf<-currDf
  }else{
    outDf<-rbind(outDf,currDf)
  }
}
write.table(outDf,file = bedName,row.names = F,col.names = F,quote = F)

#Run depth calculations
if(length(depthVec)>0){
  library(parallel)
  currCode <- paste0(
    "cd ~/../SharedUser/FastNeutron/alignments/; ",
    "samtools depth -b ",bedName," ", bamVec," > ",
    depthVec
  )
  nCores <- max(c(1,min(c(length(currCode),ceiling(detectCores()/2)))))
  cl <- makeCluster(nCores)
  parSapply(cl = cl,X = currCode,FUN = function(x){cat("\n",x,"\n");out <- system(x,wait = T);cat("\n",x," done!\n"); return(c(x,out))})
  stopCluster(cl)
}

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