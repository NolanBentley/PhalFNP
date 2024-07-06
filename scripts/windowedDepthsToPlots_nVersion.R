#### Setup environment ####
##### Variables ####
wd <- "~/Experiments/PhalFNP/"
aggFile <- "data_ignored/secondary/windowedNs5x.csv"
aggName <- "5x5Kbx40Kb"
outDir  <- "./depthImages/"

##### Load data ####
setwd(wd)
source("./scripts/interactiveAggPlot.R")
library(ggplot2)
library(hexbin)
library(plotly)
library(htmlwidgets)
aggDf <- read.csv(aggFile)

##### Format id ####
if(is.null(aggDf$id_orig)){aggDf$id_orig <- aggDf$id}
aggDf$idPre <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\1",aggDf$id)
aggDf$idNum <- as.numeric(gsub("(.*FIL..)_(\\d+)_(.*)","\\2",aggDf$id))
aggDf$idSuf <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\3",aggDf$id)
aggDf$idNumBuffered <- gsub(" ","0",format(aggDf$idNum))
aggDf$id    <- paste0(aggDf$idPre,"_",aggDf$idNumBuffered,"_",aggDf$idSuf)
aggDf$idNum_bin <- gsub(" ","0",paste0("Lines",format(floor(aggDf$idNum/100)*100),"-",format((floor(aggDf$idNum/100)+1)*100-1)))

##### Format chr ####
if(is.null(aggDf$chr_orig)){aggDf$chr_orig <- aggDf$chr}
aggDf$chrBuffered <- gsub(" ","0",gsub(".*Chr","Chr",format(aggDf$chr_orig,justify = "right")))
aggDf$chrLogic <- grepl("^(C|c)hr",aggDf$chr_orig)
aggDf$chr[!aggDf$chrLogic] <- "scaffold"
aggDf$id_chr <- paste0(aggDf$id,"_",aggDf$chr)
aggDf$id_seq <- paste0(aggDf$id,"_",aggDf$chrBuffered)

##### Find subset with divergent patterns#####
aggDens <- density(aggDf$median[aggDf$chrLogic],bw = 0.03)
aggDens$haploidLogic    <- aggDens$x< 1.5
aggDens$diploidLogic    <- aggDens$x> 1.5 & aggDens$x< 2.5
aggDens$triploidLogic   <- aggDens$x> 2.5 & aggDens$x< 3.5
aggDens$tetraploidLogic <- aggDens$x> 3.5 & aggDens$x< 4.5
aggDens$pentaploidLogic <- aggDens$x> 4.5 & aggDens$x< 5.5
aggDens$hexaploidLogic  <- aggDens$x> 5.5 & aggDens$x< 6.5

peaks <- c(
  haploid   = aggDens$x[aggDens$haploidLogic   ][which.max(aggDens$y[aggDens$haploidLogic   ])],
  diploid   = aggDens$x[aggDens$diploidLogic   ][which.max(aggDens$y[aggDens$diploidLogic   ])],
  triploid  = aggDens$x[aggDens$triploidLogic  ][which.max(aggDens$y[aggDens$triploidLogic  ])],
  tetraploid= aggDens$x[aggDens$tetraploidLogic][which.max(aggDens$y[aggDens$tetraploidLogic])],
  pentaploid= aggDens$x[aggDens$pentaploidLogic][which.max(aggDens$y[aggDens$pentaploidLogic])],
  hexaploid = aggDens$x[aggDens$hexaploidLogic ][which.max(aggDens$y[aggDens$hexaploidLogic ])]
)

aggDens$HapToDipLogic<- aggDens$x>peaks["haploid"]&aggDens$x<peaks["diploid"]
aggDens$DipToTriLogic<- aggDens$x>peaks["diploid"]&aggDens$x<peaks["triploid"]
cutoffs <- c(
  aggDens$x[aggDens$HapToDipLogic][which.min(aggDens$y[aggDens$HapToDipLogic])],
  aggDens$x[aggDens$DipToTriLogic][which.min(aggDens$y[aggDens$DipToTriLogic])]
)
png("data/densityPlot_ploidy.png",width = 8,height = 7,units = "in",res = 600)
plot(aggDens,ylim=c(0,aggDens$y[aggDens$x%in%peaks["haploid"]]*1.2))
abline(v=cutoffs,col="red",lwd=0.5)
abline(v = peaks,col="blue",lwd=1,lty=2)
lines(aggDens,lwd=2)
text(y = c(aggDens$y[aggDens$x%in%peaks["haploid"]]*1.15,aggDens$y[aggDens$x%in%peaks["haploid"]]*1.1),x = peaks,
     labels = paste0(names(peaks),"\n(",round(peaks,3),")"),
     col="navyblue",lwd=2,lty=2)
dev.off()

#Add results into aggDf
aggDf$hasDivergentMedian <- aggDf$median<=min(cutoffs)|aggDf$median>=max(cutoffs)
divPerLG <- table(aggDf$id_seq[aggDf$hasDivergentMedian])
aggDf$id_chr_inDiv <- aggDf$id_seq%in%names(divPerLG)

mean(aggDf$id_chr_inDiv)

#### Format descriptors ##### 
##### Create file names for plots ####
aggDf$scaffoldOrChr       <- gsub("Chr.*","Chr",aggDf$chr) 
aggDf$multiLineMultiChr   <- paste0(outDir,"multiChr/LSVPlot_", aggName,"_",aggDf$scaffoldOrChr,".html")
aggDf$multiLineSingleChr  <- paste0(outDir,"singleChr/LSVPlot_",aggName,"_",aggDf$chr,          ".html")
aggDf$singleLineMultiChr  <- paste0(outDir,"multiChr/",aggDf$idNum_bi,"/LSVPlot_",aggName,"_",aggDf$id,".html")
aggDf$singleLineSingleChr <- paste0(outDir,"singleChr/singleLine_",aggDf$chr,"/",aggDf$idNum_bi,"/LSVPlot_",aggName,"_",aggDf$chr,"_",aggDf$id,".html")

#### Plots #####
##### Add values for plots ####
aggDf$IntervalMidpoint_Mbp <- (aggDf$max+aggDf$min)/2000000
aggDf$IntervalRange_bp <- aggDf$max-aggDf$min+1
aggDf <- aggDf[order(aggDf$chr_orig,aggDf$chr,aggDf$idNumBuffered,aggDf$min),]

### testing plots
#testDf <- aggDf[aggDf$id%in%c("phal_FIL20_020_H_M2_1",aggDf$id[which.max(aggDf$median)]),]
#aggPlotFun(plottedDf = testDf,fileVec = testDf$singleLineSingleChr)
#aggPlotFun(plottedDf = testDf,fileVec = testDf$multiLineSingleChr)
#aggPlotFun(plottedDf = testDf,fileVec = testDf$multiLineMultiChr)
#aggPlotFun(plottedDf = testDf,fileVec = testDf$singleLineMultiChr)

##### Load gene data ####
geneDf <- read.delim("data_ignored/primary/assembly/Phallii_590_v3.2.gene.gff3",header = F,skip = 3)
colnames(geneDf)<-c("seqid","source","type","start","end","score","strand","phase","attributes")
geneDf<-geneDf[geneDf$type=="gene",]
geneDf$seqid<-gsub("scaffold_","",geneDf$seqid)
unique(geneDf$seqid)
geneDf$x   <- geneDf$start
geneDf$xend<- geneDf$end
geneDf$x   [geneDf$strand=="-"]<- geneDf$end  [geneDf$strand=="-"]
geneDf$xend[geneDf$strand=="-"]<- geneDf$start[geneDf$strand=="-"]
geneDf$overlaps1 <- c(F,geneDf$end[1:(length(geneDf$end)-1)]>geneDf$start[2:(length(geneDf$start))])
geneDf$y <- cumsum(geneDf$overlaps1)
geneDf$y_local <- zoo::rollmean(geneDf$y,k = 7,fill = "extend")
geneDf$yOffset <- (geneDf$y-geneDf$y_local)/10
hist(geneDf$yOffset,100)

##### Everything plot ####
aggPlotFun(aggDf,aggDf$multiLineMultiChr,geneDf)

### Multi-line single-chromosome plots ###
aggPlotFun(aggDf,aggDf$multiLineSingleChr,geneDf)

### Single-line multi-chromosome plots ###
aggPlotFun(aggDf,aggDf$singleLineMultiChr)

### Single-line single-chromosome plots ###
#aggPlotFun(aggDf,aggDf$singleLineSingleChr)

#### Save aggDf ####
write.csv(aggDf,"data_ignored/secondary/plottedAggDf_n.csv")

#### Making summary values ####
aggDf_presub<- aggDf_presub[order(aggDf_presub$id,aggDf_presub$chr,aggDf_presub$min),]
aggDf_presub$div <- aggDf_presub$median<min(cutoffs)|aggDf_presub$median>max(cutoffs)

aggDf_presub$divMinus2 <- c(F,F,aggDf_presub$div[1:(nrow(aggDf_presub)-2)])
aggDf_presub$divMinus1 <- c(F,aggDf_presub$div[1:(nrow(aggDf_presub)-1)])
aggDf_presub$divPlus1 <- c(aggDf_presub$div[2:(nrow(aggDf_presub))],F)
aggDf_presub$divPlus2 <- c(aggDf_presub$div[3:(nrow(aggDf_presub))],F,F)

aggDf_presub$Win3Sum <- aggDf_presub$divMinus1 + aggDf_presub$div +  aggDf_presub$divPlus1
aggDf_presub$Win5Sum <- aggDf_presub$divMinus1 + aggDf_presub$div +  aggDf_presub$divPlus1 + aggDf_presub$divMinus2 + aggDf_presub$divPlus2

win3Max <- aggregate(aggDf_presub$Win3Sum,by=list(aggDf_presub$id_seq),max)
paste0("Of Id+Seq combos: ",sum(win3Max$x>=3)," (",round(mean(win3Max$x>=3)*100,2),"% of ",nrow(win3Max)," combos)")
aggDf_noScaffold <- aggDf_presub[aggDf_presub$chrLogic,]
win3Max <- aggregate(aggDf_noScaffold$Win3Sum,by=list(aggDf_noScaffold$id_seq),max)
paste0("Of Id+chr combos (no scaffolds): ",sum(win3Max$x>=3)," (",round(mean(win3Max$x>=3)*100,2),"% of ",nrow(win3Max)," combos)")
win3Max <- aggregate(aggDf_noScaffold$Win3Sum,by=list(aggDf_noScaffold$id),max)
paste0("Of Ids: ",sum(win3Max$x>=3)," (",round(mean(win3Max$x>=3)*100,2),"% of ",nrow(win3Max)," combos)")

win5Max <- aggregate(aggDf_presub$Win5Sum,by=list(aggDf_presub$id_seq),max)
paste0("Of Id+Seq combos: ",sum(win5Max$x>=5)," (",round(mean(win5Max$x>=5)*100,2),"% of ",nrow(win5Max)," combos)")
aggDf_noScaffold <- aggDf_presub[aggDf_presub$chrLogic,]
win5Max <- aggregate(aggDf_noScaffold$Win5Sum,by=list(aggDf_noScaffold$id_seq),max)
paste0("Of Id+chr combos (no scaffolds): ",sum(win5Max$x>=5)," (",round(mean(win5Max$x>=5)*100,2),"% of ",nrow(win5Max)," combos)")
win5Max <- aggregate(aggDf_noScaffold$Win5Sum,by=list(aggDf_noScaffold$id),max)
paste0("Of Ids: ",sum(win5Max$x>=5)," (",round(mean(win5Max$x>=5)*100,2),"% of ",nrow(win5Max)," combos)")











