#Setup
wd <- "~/Experiments/PhalFNP/"; setwd(wd)

#Create read depth data
bamVecFull <- list.files("data_ignored/primary/bam/",pattern = "bam$",full.names = T)
bamVecSplit <- split(1:length(bamVecFull), cut(seq_along(1:length(bamVecFull)), breaks = 10, labels = FALSE))
i<-1
for(i in 1:length(bamVecSplit)){
  bamVec <- bamVecFull[bamVecSplit[[i]]]
  currOutFile <- "data_ignored/secondary/readCount.txt"
  currScriptFile <- paste0("scripts/readCounts/readCounting",i,".sh")
  currCmd <- c(
    paste0("echo '#Read counts of primary mapped reads' > ",currOutFile),
    as.vector(rbind(paste0("echo '",bamVec,"' >> ",currOutFile),
                    paste0("samtools view -c -F 260 ",bamVec," >> ",currOutFile))
    ))
  write(currCmd,file = currScriptFile)
  system(paste0("bash ",currScriptFile),wait = F)
  print(i)
}