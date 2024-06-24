#Setup variables
wd       <- "~/Experiments/PhalFNP"
depthDir <- "./data_ignored/secondary/depths"
outDir   <- "./data_ignored/secondary"
prelimFile <- paste0("prelimAggregateDepths_",format(Sys.time(), "%Y%m%d_%H%M"),".csv")
totalSteps <- 5
nCores     <- min(c(10,parallel::detectCores()/2))

#Setup the environment
setwd(wd)
library(parallel)

## Load functions
# Calculate the mean depth per window (byLen)
source('./scripts/aggFun1.R')

#Get data
depthFileList <- lapply(1:totalSteps,function(x){list.files(depthDir,pattern = paste0(x,"of",totalSteps,"int.depth\\.out$"),full.names = T)})
cl <- makeCluster(nCores)
j<-1
for(j in 1:length(depthFileList)){
  if(j==1){aggDf <- NULL}
  byLen <- NULL #Triggers aggFun1 to recalculate each step
  cat("\n\n\n",j,"\n")
  currDepthFiles <- depthFileList[[j]]
  #currDepthFiles<-currDepthFiles[1:40]
  parOut   <- parallel::parSapply(cl,currDepthFiles,aggFun2,simplify = F)
  parOutDf <- do.call("rbind",parOut)
  #View(parOutDf)
  aggDf <- rbind(aggDf,parOut)
}
write.csv(aggDf,file = file.path(outDir,prelimFile))
stopCluster(cl)

