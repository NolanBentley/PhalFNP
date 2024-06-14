#Setup variables
bamVec <-list.files("~/../SharedUser/FastNeutron/alignments/",pattern = "bam$",full.names = T)
byLen <- 200000
startPos <- 50
regionWidth <- 5000

initDf <- data.frame(
  chr=c('Chr01','Chr02','Chr03','Chr04','Chr05','Chr06','Chr07','Chr08','Chr09'),
  max=c('59908650','57868950','64542600','52912950','59521150','45727350','48048900','44948900','73933400')
)

#Create filenames
byLenText <- paste0(gsub("\\+","",format(byLen)),"x",regionWidth)
bedName <- paste0("~/../SharedUser/FastNeutron/alignments/regions_",byLenText,".bed")
depthVec <- paste0(gsub("\\.bam$","",bamVec),"_",byLenText,"int.depth.out")


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
library(parallel)
currCode <- paste0(
  "cd ~/../SharedUser/FastNeutron/alignments/; ",
  "samtools depth -b ",bedName," ", bamVec," > ",
  depthVec
)
nCores <- max(c(1,min(c(length(currCode),ceiling(detectCores()/2)))))
cl <- makeCluster(nCores)
applyOut <- parSapply(cl = cl,X = currCode,FUN = function(x){cat("\n",x,"\n");out <- system(x,intern = T);cat("\n",x," done!\n"); return(c(x,out))})
stopCluster(cl)

#Visualize depths
library(data.table)
outLs <- list()
for(i in 1:length(depthVec)){
  outLs[[i]] <- fread(depthVec[i])
  outLs[[i]]$marker   <- paste0(outLs[[i]]$V1,"_",outLs[[i]]$V2)
  outLs[[i]]$interval <- floor(outLs[[i]]$V2/byLen)
  outLs[[i]]$within   <- outLs[[i]]$V2%%byLen
}

lsAgg <-lapply(outLs,function(x){aggregate(x$V3,by=list(x$V1,x$interval),mean)})
for(i in 1:length(lsAgg)){
  if(i==1){aggDf <- NULL}
  aggDf <- rbind(aggDf,data.frame(
    ind = basename(depthVec[i]),
    chr = lsAgg[[i]]$Group.1,
    int = lsAgg[[i]]$Group.2,
    avg =lsAgg[[i]]$x
  ))
}
aggDf$region <- paste0(aggDf$chr,"_",aggDf$int)


aggMedian <- aggregate(aggDf$avg,by=list(aggDf$ind),median)
aggMad    <- aggregate(aggDf$avg,by=list(aggDf$ind),mad   )

aggDf$sampleMedian<-aggMedian$x[match(aggDf$ind,aggMedian$Group.1)]
aggDf$sampleMad   <-aggMad$x[match(aggDf$ind,aggMad$Group.1)]

aggDf$avg_samNorm <- (aggDf$avg-aggDf$sampleMedian)/aggDf$sampleMad

regMedian <- aggregate(aggDf$avg_samNorm,by=list(aggDf$region),median)
regMad    <- aggregate(aggDf$avg_samNorm,by=list(aggDf$region),mad   )
aggDf$regMedian<-regMedian$x[match(aggDf$region,regMedian$Group.1)]
aggDf$regMad   <-regMad$x[match(aggDf$region,regMad$Group.1)]

aggDf$avg_FullNorm <- (aggDf$avg_samNorm-aggDf$regMedian)


library(ggplot2)
ggplot(aggDf,aes(int,avg_FullNorm,color=ind,group=ind))+
  geom_line()+
  facet_wrap(~chr,ncol = 1)
 
