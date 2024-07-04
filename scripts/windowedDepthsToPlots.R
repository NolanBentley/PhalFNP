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
##### Format id ####
if(is.null(aggDf$id_orig)){aggDf$id_orig <- aggDf$id}
aggDf$idPre <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\1",aggDf$id)
aggDf$idNum <- as.numeric(gsub("(.*FIL..)_(\\d+)_(.*)","\\2",aggDf$id))
aggDf$idSuf <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\3",aggDf$id)
aggDf$idNumBuffered <- gsub(" ","0",format(aggDf$idNum))
aggDf$id    <- paste0(aggDf$idPre,"_",aggDf$idNumBuffered,"_",aggDf$idSuf)

##### Format chr ####
if(is.null(aggDf$chr_orig)){aggDf$chr_orig <- aggDf$chr}
aggDf$chrBuffered <- gsub(" ","0",gsub(".*Chr","Chr",format(aggDf$chr_orig,justify = "right")))
aggDf$chrLogic <- grepl("^(C|c)hr",aggDf$chr_orig)
aggDf$chr[!aggDf$chrLogic] <- "scaffold"
aggDf$id_chr <- paste0(aggDf$id,"_",aggDf$chr)
aggDf$id_seq <- paste0(aggDf$id,"_",aggDf$chrBuffered)

##### Find subset with divergent patterns#####
aggDens <- density(aggDf$median[aggDf$chrLogic],adjust = 0.1)
plot(aggDens,ylim=c(0,0.06))
cuttoffs <- c(-4.5,3.5)
abline(v=cuttoffs,col="red")

aggDf$hasDivergentMedian <- aggDf$median<=min(cuttoffs)|aggDf$median>=max(cuttoffs)

divPerLG <- table(aggDf$id_seq[aggDf$hasDivergentMedian])
aggDf$id_chr_inDiv <- aggDf$id_seq%in%names(divPerLG)

mean(aggDf$id_chr_inDiv)

##### Subset affDf based on having at least 1 divergent median per LG ####
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
aggDf$idNum_bin <- gsub(" ","0",paste0("Lines",
  format(floor((aggDf$idNum)/100)*100),"-",format((floor(((aggDf$idNum)-1)/100)+1)*100-1)
))

##### Create file names for plots ####
aggDf$scaffoldOrChr       <- gsub("Chr.*","Chr",aggDf$chr) 
aggDf$multiLineMultiChr   <- paste0(outDir,"multiChr/LSVPlot_", aggName,"_",aggDf$scaffoldOrChr,".html")
aggDf$multiLineSingleChr  <- paste0(outDir,"singleChr/LSVPlot_",aggName,"_",aggDf$chr,          ".html")
aggDf$singleLineMultiChr  <- paste0(outDir,"multiChr/",aggDf$idNum_bi,"/LSVPlot_",aggName,"_",aggDf$id,".html")
aggDf$singleLineSingleChr <- paste0(outDir,"singleChr/singleLine_",aggDf$chr,"/",aggDf$idNum_bi,"/LSVPlot_",aggName,"_",aggDf$chr,"_",aggDf$id,".html")

#### Build index ####
index <- c('<!DOCTYPE html>','<html>','<head>','<title>My GitHub Pages Site</title>','</head>','<body>','<h1>Hello world!</h1>','<p>Welcome to my GitHub Pages site!</p>','</body>','</html>')
index_endBody <- which(trimws(index)=="</body>")
index_added <- c(
  "<br><br><b><u>multiLineMultiChr</b></u><ul>",
  sort(unique(aggDf$multiLineMultiChr)),
  "</ul><br><br><b><u>multiLineSingleChr</b></u><ul>",
  sort(unique(aggDf$multiLineSingleChr)),
  "</ul><br><br><b><u>Sample Bins</b></u><ul>",
  sort(unique(aggDf$idNum_bin)),"/ul"
)
htmlPos <-grep("html",index_added)
index_mod <- index_added
index_mod[htmlPos] <- paste0('<li><a href="',
       gsub("^\\.","/PhalFNP",index_added[htmlPos]),
       '">',basename(index_added[htmlPos]),
       "</a></li>"
)
index_built <- c(
  index[1:(index_endBody-1)],
  index_mod,
  index[index_endBody:length(index)]
)
write(index_built,"index.html")

#### Plots #####
##### Add values for plots ####
aggDf$IntervalMidpoint_Mbp <- (aggDf$max+aggDf$min)/2000000
aggDf$IntervalRange_bp <- aggDf$max-aggDf$min+1
aggDf <- aggDf[order(aggDf$chr_orig,aggDf$chr,aggDf$idNumBuffered,aggDf$min),]

##### Everything plot ####
aggPlotFun(aggDf,aggDf$multiLineMultiChr)

### Multi-line single-chromosome plots ###
aggPlotFun(aggDf,aggDf$multiLineSingleChr)

### Single-line multi-chromosome plots ###
aggPlotFun(aggDf,aggDf$singleLineMultiChr)

### Single-line single-chromosome plots ###
aggPlotFun(aggDf,aggDf$singleLineSingleChr)

### testing plots
'
testDf <- aggDf[aggDf$id%in%c("phal_FIL20_020_H_M2_1",aggDf$id[which.max(aggDf$median)]),]
aggPlotFun(plottedDf = testDf,fileVec = testDf$singleLineSingleChr)
aggPlotFun(plottedDf = testDf,fileVec = testDf$multiLineSingleChr)
aggPlotFun(plottedDf = testDf,fileVec = testDf$multiLineMultiChr)
aggPlotFun(plottedDf = testDf,fileVec = testDf$singleLineMultiChr)
'

