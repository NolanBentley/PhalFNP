#Setup environment
wd <- "~/Experiments/PhalFNP";setwd(wd)
kVal <- 7
divCutoffs <- c(1.5,2.5)
outDir  <- "./depthImages"
unzippedCovPath <- "./data_ignored/secondary/cnv/unzipped"
winName <- "7x10kb"
badLocusSpacing <- 2500

#Load functions
source("./scripts/functions/prepareDepths.R")
source("./scripts/functions/interactiveAggPlot.R")

#Load libraries
library(ggplot2)
library(hexbin)
library(ggh4x)
library(plotly)
library(htmlwidgets)
library(zoo)

#Find files
covFiles <- list.files(path = unzippedCovPath,full.names = T)
idName   <- gsub("__.*","",gsub("\\.cov$","",basename(covFiles)))
fullPlotFiles <- file.path("data_ignored/simplePlots",paste0(idName,"_Chr00.png"))

#Determine id buffer width
id10Power <-nchar(max(as.numeric(gsub("_.*","",gsub(".*FIL.._","",idName)))))-1

#Load gene data
if(!file.exists("data/geneDf.csv")){source("scripts/functions/generateGeneDf.R")}
geneDf <- read.csv("data/geneDf.csv")

#First loop to accumulate strange loci for downstream filtering
currFileInd<-which(idName=="phal_FIL20_47_F_M2_1")
for(currFileInd in 1:length(covFiles)){
  #Load data
  indDepth <- prepareDepths(i = currFileInd,currFile = covFiles[currFileInd],currId = idName[currFileInd],k = kVal)
  indDepth <- cbind(formatId(indDepth$id,id10Power),indDepth)
  
  #Determine unusual values 
  indDepth$divLogic    <- indDepth$CN      >=max(divCutoffs)|indDepth$CN      <=min(divCutoffs)
  indDepth$winDivLogic <- (indDepth$winMedCN>=max(divCutoffs)|indDepth$winMedCN<=min(divCutoffs))&indDepth$winLogic
  
  #Add divergent rows to summary table 
  if(currFileInd==1){compDV <- NULL}
  compDV <- rbind(compDV,indDepth[indDepth$divLogic|indDepth$winDivLogic,])
  
  print(currFileInd)
}

#Analyse divergent loci for overepresented windows
i<-1
compDV$overlaps <- NA
compMat <- as.matrix(compDV[,c("i","chrN","start","end","mid","overlaps")])
uniIs <- unique(compDV$i)
i<-1
for(i in 1:length(uniIs)){
  iLogic <- compMat[,"i"]==uniIs[i]
  currIsI  <- compMat[ iLogic,]
  currNotI <- compMat[!iLogic,]
  for(j in 1:nrow(currIsI)){
    currMat  <- matrix(currIsI[j,c("chrN","mid")],nrow=nrow(currNotI),ncol=2,byrow = T)
    currIsI[j,"overlaps"]<-sum(
        currNotI[,2]==currMat[,1]&
        currNotI[,3]<=currMat[,2]&
        currNotI[,4]>=currMat[,2]
    )
  }
  compMat[ iLogic,] <- currIsI
  print(i)
}
overlapCutoff<-0.05*max(compMat[,"overlaps"],na.rm = T)
hist(compMat[,"overlaps"],1000)
abline(v = overlapCutoff,col="red")

#Save
if(!file.exists("data_ignored/secondary/compDV.csv")){
  write.csv(compDV,"data_ignored/secondary/compDV.csv",row.names = F)
}

badLoci <- as.data.frame(compMat[which(compMat[,"overlaps"]>overlapCutoff),])
badLoci <- badLoci[order(badLoci$chrN,badLoci$start,badLoci$end),]
ith <- 1
whileT <- T
while(whileT){
  currLocus <- badLoci[ith,]
  currLogic <- 
    badLoci$chrN==currLocus$chrN&
    badLoci$start-badLocusSpacing <= currLocus$mid &
    badLoci$end + badLocusSpacing >= currLocus$mid
  currLogic[(1:nrow(badLoci))<=ith]<-F
  badLoci <- badLoci[!currLogic,]
  if(ith==nrow(badLoci)){
    whileT <- F
  }
  print(paste0(ith," (",nrow(badLoci),", -",sum(currLogic),")"))
  ith <- ith+1
}
write.csv(badLoci,"data/highlyDivergentCNLoci.csv")


#Loop over files WITH FILTERING (not implemented yet)
currFileInd<-which(idName=="phal_FIL20_47_F_M2_1")
for(currFileInd in 1:length(covFiles)){
  #Load data
  indDepth <- prepareDepths(i = currFileInd,currFile = covFiles[currFileInd],currId = idName[currFileInd],k = kVal,
                            lociToRemove = badLoci,filterOnBadLoci = T)
  indDepth <- cbind(formatId(indDepth$id,id10Power),indDepth)
  
  #For checking subsetting indDepth[unique(as.vector(mapply(`:`,(which(diff(indDepth$order)!=1)-10),(which(diff(indDepth$order)!=1)+10)))),]
  
  # Create file names for plots
  indDepth$scaffoldOrChr       <- c("scaffold","Chr")[grepl("^Chr",indDepth$chr)+1] 
  indDepth$multiLineMultiChr   <- file.path(outDir,"multiChr" ,paste0("LSVPlot_", winName,"_",indDepth$scaffoldOrChr,".html"))
  indDepth$multiLineSingleChr  <- file.path(outDir,"singleChr",paste0("LSVPlot_", winName,"_",indDepth$chr,          ".html"))
  indDepth$singleLineMultiChr  <- file.path(outDir,"multiChr",indDepth$idNum_bi,paste0("LSVPlot_",winName,"_",indDepth$id_form,".html"))
  indDepth$singleLineSingleChr <- file.path(
    outDir,"singleChr",paste0("singleLine_",indDepth$chr),indDepth$idNum_bi,
    paste0("LSVPlot_",winName,"_",indDepth$chr,"_",indDepth$id_form,".html")
  )
  
  #Determine unusual values 
  indDepth$divLogic    <- indDepth$CN      >=max(divCutoffs)|indDepth$CN      <=min(divCutoffs)
  indDepth$winDivLogic <- (indDepth$winMedCN>=max(divCutoffs)|indDepth$winMedCN<=min(divCutoffs))&indDepth$winLogic
  
  #Add divergent rows to summary table 
  if(currFileInd==1){compDV2 <- NULL}
  compDV2 <- rbind(compDV2,indDepth[indDepth$divLogic|indDepth$winDivLogic,])
  
  #Make summary plots for each line
  plottedDf <- NULL
  plottedDf <- prepareCNVForAggPlot(indDepth,winName)
  if(nrow(plottedDf)>0){aggPlotFun(plottedDf,plottedDf$singleLineMultiChr,geneDf,plotIfNoDiv = T)}
  
  #Make single line and single chromosome plots only at divergent chromosomes
  plottedDf <- plottedDf[plottedDf$chr%in%unique(plottedDf$chr[plottedDf$winDivLogic]),]
  if(nrow(plottedDf)>0){aggPlotFun(plottedDf,plottedDf$singleLineSingleChr,geneDf)}
  
  #Print progress
  print(round(currFileInd/length(covFiles),4))
}

#Due to debuggining, check for repeats
compDV2$coded<-paste0(compDV2$id_form,"__",compDV2$chr,"_",compDV2$start,"..",compDV2$end)
compDV2 <- compDV2[!duplicated(compDV2$coded),]

#Make multi line plots
plottedDf <- NULL
plottedDf <- prepareCNVForAggPlot(compDV2,winName)
plottedDf <- plottedDf[plottedDf$chr%in%unique(plottedDf$chr[plottedDf$winDivLogic]),]
plottedDf <- plottedDf[order(plottedDf$chrN,plottedDf$start,plottedDf$id_form,plottedDf$end,plottedDf$CN),]
if(nrow(plottedDf)>0){aggPlotFun(plottedDf,plottedDf$multiLineMultiChr,geneDf)}
if(nrow(plottedDf)>0){aggPlotFun(plottedDf,plottedDf$multiLineSingleChr,geneDf)}

#Identify consecutive intervals
compDV2$chrNum <-as.numeric(gsub("Chr|scaffold_","",compDV2$chr))
baseLen <- nchar(max(compDV2$order))
compDV2$chrStag    <- 10^(baseLen+2)*compDV2$chrNum
compDV2$sampleStag <- 10^(baseLen+nchar(max(compDV2$chrNum))+4)*compDV2$i
compDV2$staggeredOrder <- compDV2$order+compDV2$chrStag+compDV2$sampleStag
compDV2$prevDiffs <- c(0,diff(compDV2$staggeredOrder))
compDV2$binExt    <- cumsum(abs(compDV2$prevDiffs)!=1)

#Examine binnning
visRows<-which(compDV2$binExt==as.numeric(names(which.max(table(compDV2$binExt)))))
compDV2[seq(min(visRows)-3,min(visRows)+3),]
compDV2[seq(max(visRows)-3,max(visRows)+3),]

#Save data
if(!file.exists("data_ignored/secondary/compDV2.csv")){
  write.csv(compDV2,"data_ignored/secondary/compDV2.csv",row.names = F)
}

#Check for k consecutive divergent intervals
xMax <- ceiling(max(compDV2$mid)/1000000)
geneDf <- read.csv("data/geneDf.csv")
for(kLoop in seq(3,25,by=2)){
  kVal <- kLoop
  compDV2$kSum <- zoo::rollsum(compDV2$stagDiffs,k = kVal,fill = T)
  compDV2$kIntLogic <- compDV2$kSum==kVal
  mat1 <- matrix(which(compDV2$kIntLogic),nrow = length(which(compDV2$kIntLogic)),ncol=kVal-1)
  mat2 <- matrix(seq(-1*(kVal-1)/2,(kVal-1)/2)[-((kVal-1)/2+1)],nrow = nrow(mat1),ncol=ncol(mat1),byrow=T)
  compDV2$kIntLogic[unique(as.vector(mat1+mat2))]<-T
  
  #Full "divergent" plot
  p1 <- ggplot(compDV2[compDV2$kIntLogic,],aes(mid/1000000,CN,color=id))+
    geom_hline(yintercept = 2,color="red")+
    geom_point()+
    theme_bw()+
    facet_wrap(~chr,ncol=3)+
    scale_y_continuous(limits = c(0,10  ),breaks = seq(0,10))+
    scale_x_continuous(limits = c(0,xMax),breaks = seq(0,xMax,by=5))+
    scale_color_discrete(guide="none")+
    labs(y=paste0("CN (minimum of ",kVal, " intervals = ",sum(compDV2$kIntLogic),")"),
         x="Position in Mb")+
    geom_ribbon(data = geneDf,mapping = aes(geneDf$x,y=NULL,color=NULL))
  ggsave(p1,file=paste0("data/fullCnvPlot_",kVal,".png"),width = 10,height = 8,dpi = 600)
  print(kVal)
}
