#### Setup environment ####
##### Variables ####
wd <- "~/Experiments/PhalFNP/"
aggFile <- "data_ignored/secondary/windowedDepths5x.csv"
aggName <- "5x5KbEvery40kb"
outDir  <- "./depthImages/"

##### Load data ####
setwd(wd)
source("./scripts/interactiveAggPlot.R")
library(ggplot2)
library(hexbin)
library(plotly)
library(htmlwidgets)
aggDf <- read.csv(aggFile)

#### Remove unhelpful records ####
##### Format chr ####
aggDf$chr <- gsub(" ","0",gsub(".*Chr","Chr",format(aggDf$chr)))
if(is.null(aggDf$chr_orig)){aggDf$chr_orig <- aggDf$chr}
aggDf$chrLogic <- grepl("^(C|c)hr",aggDf$chr_orig)
aggDf$chr[!aggDf$chrLogic] <- "scaffold"
aggDf$id_chr <- paste0(aggDf$id,"_",aggDf$chr)
aggDf$id_seq <- paste0(aggDf$id,"_",aggDf$chr_orig)

##### Find subset with divergent patterns#####
aggDens <- density(aggDf$median[aggDf$chrLogic],adjust = 0.1)
plot(aggDens,ylim=c(0,0.06))
cuttoffs <- c(-4.5,3.5)
abline(v=cuttoffs,col="red")

aggDf$hasDivergentMedian <- aggDf$median<=min(cuttoffs)|aggDf$median>=max(cuttoffs)

divPerLG <- table(aggDf$id_seq[aggDf$hasDivergentMedian])
aggDf$id_chr_inDiv <- aggDf$id_seq%in%names(divPerLG)

mean(aggDf$id_chr_inDiv)

if(!exists("aggDf_presub")){aggDf_presub <- aggDf}
aggDf <- aggDf[aggDf$id_chr_inDiv,]

#### Format descriptors ##### 
##### Edit line and family number ####
# To do: Replace id gsubs with split based methods
# To do: determine origin of scaffold suffixed id's 
'
idList <- strsplit(aggDf$id,split = "_")
aggDf$idSplitLen <- unlist(lapply(idList,length))
uniSplitLen <- unique(aggDf$idSplitLen)
for(i in 1:length(uniSplitLen)){
  print(head(aggDf[aggDf$idSplitLen==uniSplitLen[i],1:5]))
}
'


aggDf$id     <- gsub("Fil","FIL",aggDf$id)
aggDf$dosage <- gsub(".*(FIL..).*","\\1",aggDf$id)
aggDf$idFam <- gsub("_M2.*","",gsub(".*FIL(20|30)_","",aggDf$id))
aggDf$idNum <- as.numeric(gsub("_.*","",aggDf$idFam))

aggDf$idNum_bin <- gsub(" ","0",paste0("Lines",
  format(floor((aggDf$idNum)/100)*100),"-",format((floor(((aggDf$idNum)-1)/100)+1)*100-1)
))

aggDf$idSuffixLogic <- grepl("_",aggDf$idFam)
aggDf$idFamBuff <- paste0(aggDf$dosage,"_",gsub(" ","0",format(aggDf$idNum)))
aggDf$idFamBuff[aggDf$idSuffixLogic] <- paste0(
  aggDf$idFamBuff[aggDf$idSuffixLogic],"_",
  gsub("^.*?_(.*)","\\1",aggDf$idFam[aggDf$idSuffixLogic])
)

aggDf$idFormatted <- paste0(aggDf$idFamBuff,gsub(".*_(M.*)","\\1",aggDf$id))

##### Create file names for plots ####
aggDf$scaffoldOrChr <- gsub("Chr.*","Chr",aggDf$chr) 
aggDf$multiLineMultiChr <- paste0(outDir,"multiChr/LSVPlot_",aggName,"_",aggDf$scaffoldOrChr,".html")
aggDf$multiLineSingleChr <- paste0(outDir,"singleChr/LSVPlot_",aggName,"_",aggDf$chr,".html")
aggDf$singleLineMultiChr <- paste0(outDir,"multiChr/",aggDf$idNum_bi,"/LSVPlot_",aggName,"_",aggDf$id,".html")
aggDf$singleLineSingleChr <- paste0(outDir,"singleChr/singleLine_",aggDf$chr,"/",aggDf$idNum_bi,"/LSVPlot_",
                                  aggName,"_",aggDf$chr,"_",aggDf$id,".html")

#### Plots #####
##### Add values for plots ####
aggDf$IntervalMidpoint_Mbp <- (aggDf$max+aggDf$min)/2000000
aggDf$IntervalRange_bp <- aggDf$max-aggDf$min+1
aggDf <- aggDf[order(aggDf$chr,aggDf$idFamBuff,aggDf$min),]

##### Everything plot ####
aggPlotFun(aggDf,aggDf$multiLineMultiChr)

### Multi-line single-chromosome plots ###
aggPlotFun(aggDf,aggDf$multiLineSingleChr)

### Single-line multi-chromosome plots ###
aggPlotFun(aggDf,aggDf$singleLineMultiChr)

### Single-line single-chromosome plots ###
aggPlotFun(aggDf,aggDf$singleLineSingleChr)

### testing plots
testDf <- aggDf[aggDf$id==aggDf$id[which.max(aggDf$median)],]
aggPlotFun(testDf,testDf$singleLineSingleChr)


