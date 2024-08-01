#compress bam files
wd <- "~/Experiments/PhalFNP/"
bamPath <- "data_ignored/primary/bam/"
refPath <- "data_ignored/primary/assembly/phal_assemblyV3.fasta"

#Setup
setwd(wd)
bamFiles <- list.files(path = bamPath,pattern = "bam$",full.names = T)

#Filter if genozip files present
bamFiles <- normalizePath(bamFiles)
bamFiles <- bamFiles[!file.exists(paste0(bamFiles,".genozip"))]

#Check reference
refPath <- normalizePath(refPath)
if(!file.exists(refPath)|length(bamFiles)<1){stop("Files not found")}

# genozip --reference ../assembly/phal_assemblyV3.fasta phal_FIL30_97_M2-2__505532_1331314.bam phal_FIL30_96_M2-2__505532_1331312.bam
cmdVec <- paste0(
  "genozip --threads 10 --reference ",
  rep(refPath,length(bamFiles))," ",
  bamFiles," && rm -f ",
  bamFiles
)

#Write script to file
write(cmdVec,"data_ignored/primary/compressGenotypes.sh")
