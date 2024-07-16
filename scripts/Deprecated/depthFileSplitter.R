#Split up file
depthFile <- "~/Experiments/FastNeutron/analysis/fullDepth_20240516.out"
analysisDir <- "splitFull"

#Create filenames
firstRow <- unlist(strsplit(readLines(depthFile,1),split = "\t"))
fileDir  <- dirname(depthFile)
fileBase <- gsub("\\.out$","",basename(depthFile))


fileNames <- file.path(fileDir,analysisDir,paste0(fileBase,"_",paste0(c("chr","pos",paste0("geno",gsub(" ","0",format(1:(length(firstRow)-2),width = 4)))))))

#Generate bash scripts
baseCommands <- paste0("cut -f ",1:length(fileNames)," ",depthFile," > ",fileNames)

#Filter filenames
baseCommandsSub <- baseCommands[!file.exists(fileNames)]

#Split based on thread count
threads <- 10
threadAssignments <- rep(1:threads,ceiling(length(baseCommandsSub)/threads))[1:length(baseCommandsSub)]

#Generate individual script files
subScriptDir   <- paste0("./splittingScripts",threads)
subScriptFiles <- file.path(subScriptDir,paste0("CuttingScript_",gsub(" ","0",format(1:threads)),".sh"))
dir.create(subScriptDir)
for(i in 1:threads){
  currBashCommands <- baseCommandsSub[threadAssignments==i]
  write(currBashCommands,subScriptFiles[i])
}

#Generate run script
runScriptCommands <- paste0("bash ",paste0(subScriptFiles," &"))
write(runScriptCommands,paste0("splittingScript_n",threads,".sh"))



