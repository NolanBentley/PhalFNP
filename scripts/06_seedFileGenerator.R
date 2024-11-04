#Setup
wd <- "~/Experiments/PhalFNP/"; setwd(wd)
coverageFile <- "data/sampleDf.csv"
variantFile  <- "https://utexas.box.com/shared/static/yht8ojfafsab3btuq6xe2chw8by9onfo.csv"
seedRawFile  <- "data_ignored/primary/seed/PHAL_MTNT_Mutant Line Tracking Database_20240223Snapshot.csv"

#Load data
samDf   <- read.csv(coverageFile)
df1     <- read.csv(variantFile )
rawSeed <- read.csv(seedRawFile)

#Checking for lost lines
lostInM1 <- c("FIL30_467_M2","FIL30_730_M2","FIL30_675_M2","FIL30_686_M2","FIL30_752_M2")
samDf$sample[!gsub("dup$|b$","",samDf$sample)%in%c(gsub("\\.","_",rawSeed$ID.2),lostInM1)]
