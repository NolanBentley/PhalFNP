#Setup variables
wd <- "~/Experiments/PhalFNP/"; setwd(wd)
refPath <- "data_ignored/primary/assembly/phal_assemblyV3.fasta"
mapPath <- "data_ignored/primary/assembly/map.fa"
bamPath <- "data_ignored/primary/bam"
outPath <- "data_ignored/secondary/cnv"
scriptFileBase <- "scripts/dellyCnv/dellyCmd"
nThreads <- 10

#Create filenames
bamVec <- list.files(bamPath,pattern = "bam$",full.names = T)
bamDf <- data.frame(file=bamVec,id=basename(bamVec))
bamDf$dosage <- as.numeric(gsub(".*_FIL(..)_.*","\\1",bamDf$id))
bamDf$idNum  <- as.numeric(gsub(".*_FIL.._(.*?)_.*","\\1",bamDf$id))
bamVec <- bamVec[order(bamDf$dosage,bamDf$idNum,bamDf$id)]
baiVec <- paste0(bamVec,".bai")
covVec <- file.path(outPath,"cov",paste0(gsub("\\.bam$","",basename(bamVec)),".cov.gz"))
bcfVec <- file.path(outPath,"bcf",paste0(gsub("\\.bam$","",basename(bamVec)),".bcf"   ))
cmdVec <- paste0(scriptFileBase,gsub(" ","0",format(1:nThreads)),".sh")

#Create paths
dir.create(file.path(outPath,"cov"),showWarnings = F)
dir.create(file.path(outPath,"bcf"),showWarnings = F)
dir.create(dirname(scriptFileBase) ,showWarnings = F)

#Create cmds
indexCmd <- paste0("samtools index ",bamVec," ",baiVec)
dellyCmd <- paste0("delly cnv -a -g ",refPath," -m ",mapPath," -c ",covVec," -o ",bcfVec," ",bamVec)

#Split into threads
threadVec <- ceiling((1:length(dellyCmd))/length(dellyCmd)*(nThreads))
table(threadVec)
for(i in 1:nThreads){
  combCmd  <- as.vector(rbind(indexCmd[threadVec==i],dellyCmd[threadVec==i]))
  write(combCmd,file = cmdVec[i])
}

## Run with:
#conda create -n delly2 delly samtools
#cd ~/Experiments/PhalFNP/; conda activate delly2; bash scripts/dellyCnv/dellyCmd01.sh

