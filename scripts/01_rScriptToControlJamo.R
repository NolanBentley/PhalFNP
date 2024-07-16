#Variables
downloadDir <- "~/Experiments/PhalFNP/data_ignored/primary/bam/"
pathToCurlScript  <- "../jamoDepthRequests2.sh" #This file contains commands that download the files
simulConn <- 10 #The maxmimum number of connections to attempt at any given time

#Setup
setwd(downloadDir)
if(!file.exists(pathToCurlScript)){
  stop("Script file not found")
}

#Get cookies
system("curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=nolanbentley@utexas.edu' --data-urlencode 'password=*vGQM4b6GeGFStc' -c cookies > /dev/null")

#Load code and list of files already downloaded
codeLines <- readLines(pathToCurlScript)
codeOgFiles <- gsub(".*-o '(.*)'$","\\1",codeLines)
ogFiles <- list.files()
ogFiles <- ogFiles[ogFiles%in%codeOgFiles]

#Detect non downloaded files based on header
if(F){
  i<-1
  for(i in 1:ceiling(length(ogFiles)/100)){
    if(i==1){headOut <- NULL}
    currPos <- ((i-1)*100+1):min(c((i)*100,length(ogFiles)))
    headOut <- c(headOut,system(paste0("head -c 7 ",paste0(ogFiles[currPos],collapse = " ")),intern = T))
    print(i)
  }
  if(!(length(ogFiles)*2==length(headOut)&&all(ogFiles==gsub("==> | <==","",headOut[(1:(length(headOut)/2))*2-1])))){
    stop("Malformed headout)")
  }
  ogFilesSub <- ogFiles[!grepl("error$",as.character(headOut[(1:length(ogFiles))*2]))]
  
  codeSubset <- codeLines[!codeOgFiles%in%ogFilesSub]
}

#Detect non donwloaded files based on precense
if(T){
  ogFilesSub <- codeOgFiles[!file.exists(codeOgFiles)]
  #ogFilesSub <- ogFilesSub[grep("\\.bam$|\\.bai$",ogFilesSub)]
  
  codeSubset <- codeLines[codeOgFiles%in%ogFilesSub]
}


if(length(codeSubset)==0){
  stop("No code to run")
}else{
  cat(paste0("\n Starting to run ",length(codeSubset)," commands\n"))
}

#Run parallel curl commands
library(parallel)
batchSize <- 1000
cl <- NULL
loopName <- strftime(as.POSIXlt(Sys.time()),format = "%Y-%m-%d_%H-%M-%S")
for(i in 1:ceiling(length(codeSubset)/batchSize)){
  #Isolate lines
  currLines <- seq((i-1)*batchSize+1,min(c(i*batchSize,length(codeSubset))),by=1)
  #Check on parallel
  currCores <- max(c(1,min(c(simulConn,detectCores()-2,length(currLines)))))
  if(length(cl)!=currCores){
    cl <- makeCluster(max(c(1,min(c(simulConn,detectCores()-2,length(currLines))))))
  }
  #Create log file
  currBatchTime <- strftime(as.POSIXlt(Sys.time()),format = "%Y-%m-%d_%H-%M-%S")
  logFile <- paste0("logFile_Loop",loopName,"_Batch",currBatchTime,"_Lines",min(currLines),"-",max(currLines),".log")
  write(c("#Code subset:",codeSubset[currLines],"#### END ####","#Output:"),file = logFile)
  #Run batch
  parSapply(cl = cl,X = codeSubset[currLines],FUN = system, wait=T)
  #Write to log
  write(file = logFile,append = T,
        x = c(strftime(as.POSIXlt(Sys.time()),format = "%Y-%m-%d_%H-%M-%S"),"Done."))
  print(i)
}
stopCluster(cl)
