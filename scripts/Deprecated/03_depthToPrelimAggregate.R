#Setup variables
wd       <- "~/Experiments/PhalFNP"
depthDir <- "./data_ignored/secondary/depths"
outDir   <- "./data_ignored/secondary"
prelimFile <- paste0("prelimAggregateDepths_",format(Sys.time(), "%Y%m%d_%H%M"),".csv")
totalSteps <- 10
nCores     <- min(c(20,floor(parallel::detectCores()/4*3)))

#Setup the environment
setwd(wd)
library(parallel)
library(data.table)

## Load functions
# Calculate the mean depth per window (byLen)
source('./scripts/functions/aggFun1.R')

#Get data
depthFileList <- lapply(1:totalSteps,function(x){list.files(depthDir,pattern = paste0(x,"of",totalSteps,"int.depth\\.out$"),full.names = T)})
write.csv(unlist(depthFileList),file = paste0("examples/03_depthFileListOf",totalSteps,".csv"))
cl <- makeCluster(nCores)
j<-1
for(j in 1:length(depthFileList)){
  if(j==1){aggDf <- NULL}
  byLen <- NULL #Triggers aggFun1 to recalculate each step
  cat("\n\n\n",j,"\n")
  currDepthFiles <- depthFileList[[j]]
  parOut   <- parallel::parSapply(cl,currDepthFiles,aggFun2,simplify = F)
  parOutDf <- do.call("rbind",parOut)
  aggDf <- rbind(aggDf,parOutDf)
}
write.csv(aggDf,file = file.path(outDir,prelimFile))
stopCluster(cl)

#Check work
chkL <- list()
chkL$uniFiles <- unique(aggDf$file)
chkL$uniFiles_sub <- uniFiles[grep("_9._",chkL$uniFiles)]
chkL$uniFiles_sM  <- match(chkL$uniFiles_sub,aggDf$file)
chkL$subDf <- aggDf[chkL$uniFiles_sM,]
chkL$subDf$step <- as.numeric(gsub("of.*","",gsub(".*_|int\\.depth\\.out$","",chkL$subDf$file)))
View(chkL$subDf)
