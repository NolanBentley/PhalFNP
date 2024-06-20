#Setup variables
wd    <-"~/Experiments/PhalFNP"
bamPath <- "data_ignored/primary/bam"
bedPath <- "data_ignored/secondary/bed"
depthPath <- "data_ignored/secondary/depths"
fastaPath <- "data_ignored/primary/assembly/phal_assemblyV3.fasta"
byLen <- 20000
buffer <- 2000
regionWidth <- 5000
steps       <- 5
nCores <- 10

#Setup packages
packageDir <- file.path(wd,"R")
packageVec <- c("ggplot2","chron","data.table")
library(Biostrings) #Installation somewhat tricky

#Use paths
setwd(wd)
dir.create(packageDir,showWarnings = F)
.libPaths(packageDir)
if(!all(packageVec%in%installed.packages())){install.packages(packageVec,lib=packageDir, repos='http://cran.us.r-project.org')}
library("chron")

#Describe genome
ref <- Biostrings::readDNAStringSet(fastaPath)
initDf <- data.frame(
  chr=names(ref),
  max=nchar(ref)
)

#Find bam files
bamVec<-file.path(getwd(),list.files(bamPath,pattern = "bam$",full.names = T))
bamVec<-bamVec[order(as.numeric(gsub("_.*","",gsub("^.._","",gsub(".*FIL","",basename(bamVec))))))]
trimmedBam <- gsub("\\.bam","",bamVec)

#Create filenames
byLenText <- paste0(buffer/1000,"Kbuffer_",byLen/1000,"Kint_",regionWidth/1000,"Kwide_",1:steps,"of",steps)
bedName   <- file.path(getwd(),bedPath,paste0("regions_",byLenText,".bed"))
fileDf <- data.frame(bam = bamVec)
for(i in 1:length(byLenText)){
  fileDf[[i+1]]<-file.path(getwd(),depthPath,paste0(basename(trimmedBam),"_",byLenText[i],"int.depth.out"))
}
fileDf <- data.frame(bam = fileDf$bam,depth=as.vector(unlist(fileDf[,2:ncol(fileDf)])))

#Remove files with the files already extant
fileDf$exist <- file.exists(fileDf$depth)
bamVec   <- fileDf$bam[!fileDf$exist]
depthVec <- fileDf$depth[!fileDf$exist]

#Summarize fbam files analyzed
write.csv(fileDf,"data/bamFilesInDepth.csv")


#Generate bed files
for(i in 1:nrow(initDf)){
  currDf <- data.frame(chr=initDf$chr[i],min=seq(buffer,initDf[i,2]-buffer,by=byLen),
                       max=seq(buffer,initDf[i,2]-buffer,by=byLen)+regionWidth)
  if(i==1){outDf<-currDf
  }else{outDf<-rbind(outDf,currDf)}
}
outDf_split<-outDf
stepVec <- rep(1:steps,ceiling(nrow(outDf_split)/steps))[1:nrow(outDf_split)]
for(i in 1:length(bedName)){
  write.table(outDf[stepVec==i,],file = bedName[i],row.names = F,col.names = F,quote = F)
}

#Run depth calculations
if(length(depthVec)>0){
  library(parallel)
  currCode <- paste0(
    "cd ~/../SharedUser/FastNeutron/alignments/; ",
    "samtools depth -b ",bedName," ", bamVec," > ",
    depthVec
  )
  cl <- makeCluster(nCores)
  parSapply(cl = cl,X = currCode,FUN = function(x){cat("\n",x,"\n");out <- system(x,wait = T);cat("\n",x," done!\n"); return(c(x,out))})
  stopCluster(cl)
}
