#Setup environment
wd <- "~/Experiments/PhalFNP/";setwd(wd)
library(ggplot2)

#Find files
covFiles <- list.files(path = "data_ignored/secondary/cnv/unzipped/",full.names = T)
idName   <- gsub("__.*","",gsub("\\.cov$","",basename(covFiles)))
fullPlotFiles <- file.path("data_ignored/simplePlots",paste0(idName,"_Chr00.png"))

#Loop over files
for(i in 1:length(covFiles)){
  if(i==1){compDV <- NULL}
  #Read and clean
  indDepth   <- read.delim(covFiles[i])
  colnames(indDepth)<- c("chr","start","end","mappable","counts","CN")
  indDepth <- indDepth[grep("Chr",indDepth$chr),]
  indDepth$i    <- i
  indDepth$file <- covFiles[i]
  indDepth$id   <- idName[i]
  indDepth$order<- 1:nrow(indDepth)
  indDepth$mid  <- (indDepth$start+indDepth$end)/2
  
  
  #Add in unusual values to divergent data.frame
  compDV <- rbind(compDV,indDepth[indDepth$CN>=2.5|indDepth$CN<=1.5,])
  
  #Plot full
  maxX <- ceiling(max(indDepth$mid/1000000))
  p1 <- ggplot(indDepth,aes(mid/1000000,CN))+
    geom_hline(yintercept = 2,color="red")+
    geom_bin_2d(binwidth=c(0.2,0.1))+theme_bw()+
    scale_fill_viridis_c(trans="log10")+
    scale_y_continuous(limits = c(0,8   ),breaks = 0:10)+
    scale_x_continuous(limits = c(0,maxX),breaks = seq(0,maxX,by=1))+
    facet_wrap(~chr,ncol = 1)+
    labs(title=gsub("\\.cov$","",idName[i]),x="Midpoint (Mb)","Copy Number")
  ggsave(plot = p1,filename = fullPlotFiles[i],height = 26,width = 15,dpi = 600)
  
  #Print progress
  cat(i)
  cat(": ")
  cat(idName[i])
  
  #Plot singles
  j<-1
  uniChr <- unique(indDepth$chr)
  pList <- list()
  singleChrFiles <- paste0(gsub("_Chr00\\.png$","",fullPlotFiles[i]),"_",uniChr,".png")
  for(j in 1:length(uniChr)){
    pList[[j]] <- ggplot(indDepth[indDepth$chr==uniChr[j],],
                         aes(mid/1000000,CN))+
      geom_hline(yintercept = 2,color="red")+
      geom_bin_2d(binwidth=c(0.2,0.1))+theme_bw()+
      scale_fill_viridis_c(trans="log10")+
     scale_y_continuous(limits = c(0,8),breaks = 0:10)+
      scale_x_continuous(limits = c(0,maxX),breaks = seq(0,maxX,by=1))+
      facet_wrap(~chr,ncol = 1)+
      labs(title=gsub("\\.cov$","",idName[i]),x="Midpoint (Mb)","Copy Number")
    ggsave(plot = pList[[j]],filename = singleChrFiles[j],height = 8.5,width = 15,dpi = 600)
    cat(".")
  }
  cat("!\n")
}

#Identify consecutive intervals
compDV$chrNum <-as.numeric(gsub("Chr|scaffold_","",compDV$chr))
baseLen <- nchar(max(compDV$order))
compDV$chrStag    <- 10^(baseLen+2)*compDV$chrNum
compDV$sampleStag <- 10^(baseLen+nchar(max(compDV$chrNum))+4)*compDV$i
compDV$staggeredOrder <- compDV$order+compDV$chrStag+compDV$sampleStag
compDV$prevDiffs <- c(0,diff(compDV$staggeredOrder))
compDV$binExt    <- cumsum(abs(compDV$prevDiffs)!=1)
head(compDV[compDV$binExt==as.numeric(names(which.max(table(compDV$binExt)))),])

#Save or load data
if(!file.exists("data_ignored/secondary/compDV.csv")){
  write.csv(compDV,"data_ignored/secondary/compDV.csv",row.names = F)
}
if(!exists("compCV")){
  compDV <- read.csv("data_ignored/secondary/compDV.csv")
}

#Check for k consecutive divergent intervals
xMax <- ceiling(max(compDV$mid)/1000000)
geneDf <- read.csv("data/geneDf.csv")
for(kLoop in seq(3,25,by=2)){
  kVal <- kLoop
  compDV$kSum <- zoo::rollsum(compDV$stagDiffs,k = kVal,fill = T)
  compDV$kIntLogic <- compDV$kSum==kVal
  mat1 <- matrix(which(compDV$kIntLogic),nrow = length(which(compDV$kIntLogic)),ncol=kVal-1)
  mat2 <- matrix(seq(-1*(kVal-1)/2,(kVal-1)/2)[-((kVal-1)/2+1)],nrow = nrow(mat1),ncol=ncol(mat1),byrow=T)
  compDV$kIntLogic[unique(as.vector(mat1+mat2))]<-T
  
  #Full "divergent" plot
  p1 <- ggplot(compDV[compDV$kIntLogic,],aes(mid/1000000,CN,color=id))+
    geom_hline(yintercept = 2,color="red")+
    geom_point()+
    theme_bw()+
    facet_wrap(~chr,ncol=3)+
    scale_y_continuous(limits = c(0,10  ),breaks = seq(0,10))+
    scale_x_continuous(limits = c(0,xMax),breaks = seq(0,xMax,by=5))+
    scale_color_discrete(guide="none")+
    labs(y=paste0("CN (minimum of ",kVal, " intervals = ",sum(compDV$kIntLogic),")"),
         x="Position in Mb")+
    geom_ribbon(data = geneDf,mapping = aes(geneDf$x,y=NULL,color=NULL))
  ggsave(p1,file=paste0("data/fullCnvPlot_",kVal,".png"),width = 10,height = 8,dpi = 600)
  print(kVal)
}
