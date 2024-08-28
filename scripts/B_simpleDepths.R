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
#install.packages(c("ggplot2","hexbin","ggh4x","plotly","htmlwidgets","zoo"))
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


