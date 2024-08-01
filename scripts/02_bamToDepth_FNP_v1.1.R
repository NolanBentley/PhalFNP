#Setup variables
wd    <-"~/Experiments/PhalFNP"
bamPath <- "data_ignored/primary/bam"
bedPath <- "data_ignored/secondary/bed"
depthPath <- "data_ignored/secondary/depths"
fastaPath <- "data_ignored/primary/assembly/phal_assemblyV3.fasta"
byLen <- 23750
buffer <- 2000
regionWidth <- 5000
steps       <- 10
nCores <- 10

#Setup packages
packageDir <- file.path(wd,"R")
packageVec <- c("ggplot2","chron","data.table")
library(Biostrings) #Installation somewhat tricky on my server

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

i<-1
for(i in 1:length(byLenText)){
  currDepthOuts <- file.path(getwd(),depthPath,paste0(basename(trimmedBam),"_",byLenText[i],"int.depth.out"))
  if(i==1){fileDf<-NULL}
  fileDf <- rbind(fileDf,data.frame(i=i,bam=bamVec,depth=currDepthOuts))
}

#Remove files with the files already extant
fileDf$exist <- file.exists(fileDf$depth)

#Generate bed files
i<-1
for(i in 1:nrow(initDf)){
  currDf     <- data.frame(chr=initDf$chr[i],min=seq(buffer,initDf[i,2]-buffer-regionWidth,by=byLen)+1)
  currDf$max <- currDf$min+regionWidth
  currDf     <- currDf[currDf$max<=(initDf[i,2]-buffer),]
  if(i==1){outDf<-currDf
  }else if(nrow(currDf)>1){
    outDf<-rbind(outDf,currDf)
  }
}

stepVec <- rep(1:steps,ceiling(nrow(outDf)/steps))[1:nrow(outDf)]
bedName   <- file.path(getwd(),bedPath,paste0("regions_",byLenText,".bed"))
i<-1
for(i in 1:length(bedName)){
  write.table(outDf[stepVec==i,],file = bedName[i],row.names = F,col.names = F,quote = F)
}

#Add in bed files
fileDf$bed <- NA
for(i in 1:length(bedName)){
  fileDf$bed[fileDf$i==i]<-bedName[i]
}

#Generate code
fileDf$code <- paste0(
  "samtools depth -b ",fileDf$bed," ", fileDf$bam," > ",fileDf$depth
)

#Summarize bam files analyzed
## Assigned core
fileDf$core <- rep(1:nCores,ceiling(nrow(fileDf)/nCores))[1:nrow(fileDf)]
write.csv(fileDf,paste0("data/bamFilesInDepth_",regionWidth/1000,"KbEvery",byLen/1000,"KbFrom",buffer/1000,"Kb.csv"))

#Run depth calculations
i<-1
dir.create("data_ignored/secondary/depthCode")
for(i in 1:nCores){
  currCode <- fileDf$code[!fileDf$exist&fileDf$core==i]
  currScript <- paste0("data_ignored/secondary/depthCode/depthCode_",i,".sh")
  write(currCode,file = currScript)
  if(i==1){outL<-list()}
  outL[[i]] <- system(paste0("bash ",currScript),wait = F)
}