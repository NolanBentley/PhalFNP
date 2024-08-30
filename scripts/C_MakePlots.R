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

#Load gene data
if(!file.exists("data/geneDf.csv")){source("scripts/functions/generateGeneDf.R")}
geneDf <- read.csv("data/geneDf.csv")

#Load bad locus data
badLoci <- read.csv("data/highlyDivergentCNLoci.csv")

#Find files
covFiles <- list.files(path = unzippedCovPath,full.names = T)

#Testing
if(exists("debugT")&&debugT){
  covFiles<- covFiles[grep("FIL30_43_",covFiles)]
  loopVec<-2:4
}

#Determine id buffer width
idName        <- gsub("__.*","",gsub("\\.cov$","",basename(covFiles)))
fullPlotFiles <- file.path("data_ignored/simplePlots",paste0(idName,"_Chr00.png"))
id10Power <-nchar(max(as.numeric(gsub("_.*","",gsub(".*FIL.._","",idName)))))-1

#Loop over files WITH FILTERING (not implemented yet)
loopVec <- 1:length(covFiles)
for(currFileInd in loopVec){
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
  indDepth$winDivLogic <- (indDepth$winMedCN>=max(divCutoffs) | indDepth$winMedCN<=min(divCutoffs)) & indDepth$winLogic
  indDepth$divLogic    <- indDepth$CN>=max(divCutoffs) | indDepth$CN<=min(divCutoffs) | indDepth$winDivLogic
  
  #Add divergent rows to summary table 
  if(currFileInd==min(loopVec)){compDV2 <- NULL}
  compDV2 <- rbind(compDV2,indDepth[indDepth$divLogic|indDepth$winDivLogic,])
  
  #Make summary plots for each line
  plottedDf <- NULL
  plottedDf <- prepareCNVForAggPlot(indDepth,winName)
  if(nrow(plottedDf)>0){
    aggPlotFun(plottedDf = plottedDf,fileVec = plottedDf$singleLineMultiChr,genes = geneDf,plotIfNoDiv = T)
  }
  
  #Make single line and single chromosome plots only at divergent chromosomes
  plottedDf <- plottedDf[plottedDf$chr%in%unique(plottedDf$chr[plottedDf$winDivLogic]),]
  if(nrow(plottedDf)>0){aggPlotFun(plottedDf,plottedDf$singleLineSingleChr,geneDf)}
  
  #Print progress
  print(round(currFileInd/length(covFiles),4))
}

#Due to debuggining, check for repeats
compDV2$coded<-paste0(compDV2$id_form,"__",compDV2$chr,"_",compDV2$start,"..",compDV2$end)
compDV2 <- compDV2[!duplicated(compDV2$coded),]

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

#Make multi line plots
plottedDf <- NULL
plottedDf <- prepareCNVForAggPlot(compDV2,winName)
plottedDf <- plottedDf[plottedDf$chr%in%unique(plottedDf$chr[plottedDf$winDivLogic]),]
plottedDf <- plottedDf[order(plottedDf$chrN,plottedDf$start,plottedDf$id_form,plottedDf$end,plottedDf$CN),]
if(nrow(plottedDf)>0){aggPlotFun(plottedDf,plottedDf$multiLineMultiChr,geneDf)}
if(nrow(plottedDf)>0){aggPlotFun(plottedDf = plottedDf,fileVec = plottedDf$multiLineSingleChr,genes = geneDf)}

#Save an image
#save.image("data_ignored/secondary/afterCplotting.rimage")
#load("data_ignored/secondary/afterCplotting.rimage")
