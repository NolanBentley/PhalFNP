setwd("~/Experiments/PhalFNP/")
setwd("data_ignored/secondary/depths/")
dir.create("renamed")

depthList    <- list.files(pattern = "out$",full.names = T)
for(i in 1:10){
  currFiles <- depthList[grep(paste0("_",i,"of"),depthList)]
  currHeaders <- sapply(currFiles,read.delim,nrows = 5,header = F,simplify = F)
  currHeaders2 <- mapply(function(x,y){x$file <- y; return(x)},currHeaders,names(currHeaders),SIMPLIFY = F)
  currDf <- data.frame(i,do.call("rbind",currHeaders2))
  if(i==1){outDf<-NULL}
  outDf <- rbind(outDf,currDf)
  print(i)
}

#Categorize
outDf$of <- round((outDf$V2-2001)/23750+1,2)
outDf$judge <- outDf$of!=outDf$i
outDf$judge_col <- c("black","red")[(outDf$judge)+1]
outDf$order <- 1:nrow(outDf)
plot(x=outDf$order,y = outDf$V2,col=outDf$judge_col)

#Determine correct name
outDf$file_calc <- paste0(gsub("_.of.*|_..of.*","",outDf$file),"_",outDf$of,"of10int.depth.out")
rownames(outDf)<-NULL
matchVec <- match(depthList,outDf$file_calc)
unmatched <- depthList[is.na(matchVec)]

#Determine names to swtich
switchDf <- outDf[matchVec[!is.na(matchVec)],]
switchDf <- switchDf[switchDf$judge,]

#Move them
## Toggled off for safety
cmdVec <- paste0("mv ",switchDf$file," ./renamed/",gsub("^\\.\\/","",switchDf$file_calc))
if(F){
  write(cmdVec,"mv.sh")
  system("bash mv.sh")
}

#Remove the duplicated wrong files
delDf <- outDf[outDf$judge&(!outDf$file%in%switchDf$file),]
delDf <- delDf[!duplicated(delDf$file),]
if(F){
  write(paste0("rm ",delDf$file),"rm.sh")
  system("bash rm.sh")
}