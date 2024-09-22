#Setup paths
wd <- "~/Experiments/PhalFNP/"; setwd(wd)
cnvPath <- "data_ignored/secondary/cnv/unzipped"
combinedCovFile <-normalizePath("data_ignored/secondary/cnv/combined.tsv")

#Combine the files with filename delimiters
cnvFiles <- list.files(cnvPath,full.names = T)
cmdVec <- paste0(
  "echo 'title:",cnvFiles,"' > ",combinedCovFile,"; ",
  "tail --lines=+2 ",cnvFiles," >> ",combinedCovFile
)
cmdVec[-1] <- gsub(" > "," >> ",cmdVec[-1])
cmdVec[1] <- gsub("tail --lines=+2","cat",cmdVec[1],fixed = T)
cmdVec[1:2]
for(i in 1:length(cmdVec)){
  system(cmdVec[i],wait = T)
  if(i%%100==1){print(i)}
}

#Load data
library(data.table)
covDf <- fread(combinedCovFile,header=T,fill = T,sep = "\t",skip = 1)
colnames(covDf) <- gsub(".*_","",colnames(covDf))
covDf <- rbind(data.frame(chr=readLines(combinedCovFile,n = 1),start=NA,end=NA,mappable=NA,counts=NA,CN=NA),covDf)

#Fill in titles
titleRows <-grep("title\\:",covDf$chr)
titleEnds <- c(titleRows[2:length(titleRows)]-1,nrow(covDf))
id <- gsub("-.$|_.$","",gsub("phal_|__.*","",basename(gsub("title\\:","",covDf$chr[titleRows]))))
if(any(duplicated(id))){stop("ID error!")}
covDf$id <- NA
for(i in 1:length(titleRows)){
  covDf$id[(titleRows[i]):(titleEnds[i])] <- id[i]
  print(i)
}
covDf <- covDf[-titleRows,]


#Add in middle position
i<-1
j<-1
uniChr <- unique(covDf$chr)
covDf$mid <- (covDf$end+covDf$start)/2
hist(covDf$CN,10000,ylim=c(0,10000),xlim=c(0,10))
uniInds <- unique(covDf$id)
for(i in 1:length(uniChr)){
  currDf <- covDf[covDf$chr==uniChr[i],]
  currStarts <- seq(1,max(currDf$end-2500),by=2500)
  currEnds   <- currStarts+2500-1
  currMids   <- round((currStarts + currEnds)/2)
  currMatCN  <- matrix(NA,nrow = length(currMids),ncol=length(uniInds),dimnames = list(paste0(uniChr[i],"_",currMids),uniInds))
  currMatMappable <- currMatCN
  currMatCounts   <- currMatCN
  for(j in 1:length(uniInds)){
    indDf <- currDf[currDf$id==uniInds[j],]
    overlaps <- sapply(currMids,function(x){which.min(x>=indDf$start&x<indDf$end)})
    overLens <- unlist(lapply(overlaps,length))
    if(any(overLens>1)){stop("Error!")}
    overlaps[overLens==0]<-NA
    currMatCN[,uniInds[j]]<-indDf$CN[unlist(overlaps)]
    currMatMappable[,uniInds[j]]<-indDf$mappable[unlist(overlaps)]
    currMatCounts[,uniInds[j]]<-indDf$counts[unlist(overlaps)]
    print(j)
  }
  if(i==1){
    outCn  <- NULL
    outMap <- NULL
    outCnt <- NULL
  }
  outCn  <- rbind(outCn,currMatCN)
  outMap <- rbind(outMap,currMatMappable)
  outCnt <- rbind(outCnt,currMatCounts)
}
save.image("data_ignored/secondary/afterCmbCnv.rimage")

rowMedians<-apply(outCn,MARGIN = 1,median,na.rm=T)
colMedians<-apply(outCn,MARGIN = 2,median,na.rm=T)
rowMad    <-apply(outCn,MARGIN = 1,mad   ,na.rm=T)
colMad    <-apply(outCn,MARGIN = 2,mad   ,na.rm=T)
rowObs    <-apply(outCn,MARGIN = 1,function(x){sum(!is.na(x))})
colObs    <-apply(outCn,MARGIN = 2,function(x){sum(!is.na(x))})
minObs    <- 0.99*ncol(outCn)
medLims   <- c(1.8,2.2)
hist(rowObs    ,10000,ylim=c(0,200));abline(v = minObs,col="red")
hist(rowMedians,100,xlim=c(1.5,4),ylim=c(0,200));abline(v=medLims,col="red")

outCn2 <- outCn[rowMedians>=medLims[1]&rowMedians<=medLims[2]&rowObs>=minObs,]

posVec <-as.numeric(gsub(".*_","",rownames(outCn2)))
chrVec <-as.numeric(gsub("_.*","",gsub("scaffold_|Chr","",rownames(outCn2))))

currLag <- 41
library(zoo)
chrCheck  <- rollmean  (chrVec                 ,k = currLag)
posCheck  <- rollmean  (posVec                 ,k = currLag)
currCheck <- rollmedian(outCn2[,"FIL30_653_M2"],k = currLag)
chrLogic  <- (chrCheck%%1)==0

library(ggplot2)
plottedDf <- data.frame(chr=chrCheck,pos=posCheck,cn=currCheck)[chrLogic&chrCheck==9,]
ggplot(plottedDf,aes(pos,cn))+
  geom_smooth()+
  facet_wrap(~chr,ncol = 3)

unique(diff(posVec,lag = currLag)[diff(chrVec,lag = currLag)==0],10000)


#Generate place to store files
dir.create(file.path(dirname(combinedCovFile),"combSingleChr"),showWarnings = F)
outFileVec <- file.path(dirname(combinedCovFile),"combSingleChr",paste0("combCov_",uniChr,".csv"))

#Write one at a time
for(i in 1:length(uniChr)){
  print(i)
  write.csv(covDf[covDf$chr==uniChr[i],],file = outFileVec[i])
  cat("!")
}

