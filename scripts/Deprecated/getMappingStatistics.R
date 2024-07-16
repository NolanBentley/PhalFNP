#Setup
wd <- "~/Experiments/PhalFNP/"
setwd(wd)
library(parallel)

#Get files
bamVec <- list.files("data_ignored/primary/bam/",pattern = "bam$",full.names = T)[1:2]

#Run analysis
cl <- makeCluster(5)
dir.create("data_ignored/secondary/flagstatOut",showWarnings = F)
codeVec <- paste0("samtools flagstat ",bamVec," > data_ignored/secondary/flagstatOut/",gsub("\\.bam","",basename(bamVec)),"_flagstat.out")
out <- parSapply(cl,codeVec,FUN = system)
stopCluster(cl)

#Save output
saveRDS(out,file = "data_ignored/secondary/flagstatOut.rds")
