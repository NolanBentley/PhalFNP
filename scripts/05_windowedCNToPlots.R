#### Setup environment ####
##### Variables ####
wd <- "~/Experiments/PhalFNP/"
winFile <- "data_ignored/secondary/windowedNs5x.csv"
outDir  <- "./depthImages/"

##### Load data ####
setwd(wd)
source("./scripts/functions/interactiveAggPlot.R")
pcks <- c("ggplot2","hexbin","plotly","htmlwidgets","ggh4x")
if(!all(pcks%in%installed.packages())){install.packages(pcks)}
library(ggplot2)
library(hexbin)
library(plotly)
library(htmlwidgets)
library(ggh4x)#install.packages("ggh4x")
winDf <- read.csv(winFile)

winName <- paste0(unique(winDf$winName),collapse = "|")
if(winName=="7x4.999kbx23.75kb(~147.499Kb)"){
  winName <- "7x5Kbx23.75Kb"
}


##### Format id ####
if(is.null(winDf$id_orig)){winDf$id_orig <- winDf$id}
winDf$idPre <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\1",winDf$id)
winDf$idNum <- as.numeric(gsub("(.*FIL..)_(\\d+)_(.*)","\\2",winDf$id))
winDf$idSuf <-            gsub("(.*FIL..)_(\\d+)_(.*)","\\3",winDf$id)
winDf$idNumBuffered <- gsub(" ","0",format(winDf$idNum))
winDf$id    <- paste0(winDf$idPre,"_",winDf$idNumBuffered,"_",winDf$idSuf)
winDf$idNum_bin <- gsub(" ","0",paste0("Lines",format(floor(winDf$idNum/100)*100),"-",format((floor(winDf$idNum/100)+1)*100-1)))

##### Format chr ####
if(is.null(winDf$chr_orig)){winDf$chr_orig <- winDf$chr}
winDf$chrBuffered <- gsub(" ","0",gsub(".*Chr","Chr",format(winDf$chr_orig,justify = "right")))
winDf$chrLogic <- grepl("^(C|c)hr",winDf$chr_orig)
winDf$chr[!winDf$chrLogic] <- "scaffold"
winDf$id_chr <- paste0(winDf$id,"_",winDf$chr)
winDf$id_seq <- paste0(winDf$id,"_",winDf$chrBuffered)

##### Find subset with divergent patterns#####
winDens <- density(winDf$median[winDf$chrLogic],bw = 0.03)
winDens$haploidLogic    <- winDens$x< 1.5
winDens$diploidLogic    <- winDens$x> 1.5 & winDens$x< 2.5
winDens$triploidLogic   <- winDens$x> 2.5 & winDens$x< 3.5
winDens$tetraploidLogic <- winDens$x> 3.5 & winDens$x< 4.5
winDens$pentaploidLogic <- winDens$x> 4.5 & winDens$x< 5.5
winDens$hexaploidLogic  <- winDens$x> 5.5 & winDens$x< 6.5

peaks <- c(
  haploid   = winDens$x[winDens$haploidLogic   ][which.max(winDens$y[winDens$haploidLogic   ])],
  diploid   = winDens$x[winDens$diploidLogic   ][which.max(winDens$y[winDens$diploidLogic   ])],
  triploid  = winDens$x[winDens$triploidLogic  ][which.max(winDens$y[winDens$triploidLogic  ])],
  tetraploid= winDens$x[winDens$tetraploidLogic][which.max(winDens$y[winDens$tetraploidLogic])],
  pentaploid= winDens$x[winDens$pentaploidLogic][which.max(winDens$y[winDens$pentaploidLogic])],
  hexaploid = winDens$x[winDens$hexaploidLogic ][which.max(winDens$y[winDens$hexaploidLogic ])]
)

winDens$HapToDipLogic<- winDens$x>peaks["haploid"]&winDens$x<peaks["diploid"]
winDens$DipToTriLogic<- winDens$x>peaks["diploid"]&winDens$x<peaks["triploid"]
cutoffs <- c(
  winDens$x[winDens$HapToDipLogic][which.min(winDens$y[winDens$HapToDipLogic])],
  winDens$x[winDens$DipToTriLogic][which.min(winDens$y[winDens$DipToTriLogic])]
)
png("data/densityPlot_ploidy.png",width = 8,height = 7,units = "in",res = 600)
plot(winDens,ylim=c(0,winDens$y[winDens$x%in%peaks["haploid"]]*1.2))
abline(v=cutoffs,col="red",lwd=0.5)
abline(v = peaks,col="blue",lwd=1,lty=2)
lines(winDens,lwd=2)
text(y = c(winDens$y[winDens$x%in%peaks["haploid"]]*1.15,winDens$y[winDens$x%in%peaks["haploid"]]*1.1),x = peaks,
     labels = paste0(names(peaks),"\n(",round(peaks,3),")"),
     col="navyblue",lwd=2,lty=2)
dev.off()

#Add results into winDf
winDf$hasDivergentMedian <- winDf$median<=min(cutoffs)|winDf$median>=max(cutoffs)
divPerLG <- table(winDf$id_seq[winDf$hasDivergentMedian])
winDf$id_chr_inDiv <- winDf$id_seq%in%names(divPerLG)

#### Format descriptors ##### 
##### Create file names for plots ####
winDf$scaffoldOrChr       <- gsub("Chr.*","Chr",winDf$chr) 
winDf$multiLineMultiChr   <- paste0(outDir,"multiChr/LSVPlot_", winName,"_",winDf$scaffoldOrChr,".html")
winDf$multiLineSingleChr  <- paste0(outDir,"singleChr/LSVPlot_",winName,"_",winDf$chr,          ".html")
winDf$singleLineMultiChr  <- paste0(outDir,"multiChr/",winDf$idNum_bi,"/LSVPlot_",winName,"_",winDf$id,".html")
winDf$singleLineSingleChr <- paste0(outDir,"singleChr/singleLine_",winDf$chr,"/",winDf$idNum_bi,"/LSVPlot_",winName,"_",c("","scaffold")[(winDf$chr=="scaffold")+1],winDf$chr_orig,"_",winDf$id,".html")

#### Plots #####
##### Add values for plots ####
winDf$IntervalMidpoint_Mbp <- (winDf$max+winDf$min)/2000000
winDf$IntervalRange_bp <- winDf$max-winDf$min+1
winDf <- winDf[order(winDf$chr_orig,winDf$chr,winDf$idNumBuffered,winDf$min),]

##### Load gene data ####
if(!file.exists("data/geneDf.csv")){source("scripts/functions/generateGeneDf.R")}
geneDf <- read.csv("data/geneDf.csv")

### testing plots
#testDf <- winDf[winDf$id%in%c("phal_FIL20_020_H_M2_1",winDf$id[which.max(winDf$median)]),]
#aggPlotFun(plottedDf = testDf,fileVec = testDf$singleLineSingleChr,geneDf)
#aggPlotFun(plottedDf = testDf,fileVec = testDf$multiLineSingleChr,geneDf)
#aggPlotFun(plottedDf = testDf,fileVec = testDf$multiLineMultiChr,geneDf)
#aggPlotFun(plottedDf = testDf,fileVec = testDf$singleLineMultiChr,geneDf)


##### Everything plot ####
aggPlotFun(winDf,winDf$multiLineMultiChr,geneDf)

### Multi-line single-chromosome plots ###
aggPlotFun(winDf,winDf$multiLineSingleChr,geneDf)

### Single-line multi-chromosome plots ###
winDf_sub  <- winDf[winDf$id %in%(winDf$id[winDf$hasDivergentMedian])&winDf$chrLogic,]
aggPlotFun(winDf_sub,winDf_sub$singleLineMultiChr,geneDf)  

### Single-line single-chromosome plots ###
winDf_sub2 <-winDf[winDf$id_seq %in%(winDf$id_seq[winDf$hasDivergentMedian]),]
aggPlotFun(winDf_sub2,winDf_sub2$singleLineSingleChr,geneDf)

#### Save winDf ####
write.csv(winDf,"data_ignored/secondary/plottedwinDf_n.csv")
write.csv(winDf[winDf$hasDivergentMedian,],"data_ignored/secondary/winDf_DivergentOnly.csv")
#winDf <- read.csv("data_ignored/secondary/plottedwinDf_n.csv")

#Investigate sd
nonDivDf <- winDf[!winDf$hasDivergentMedian,]
win_nonDivDf <- winregate(nonDivDf$mean,by=list(nonDivDf$id),mad)

##Add in peak values
sampleDf <- read.csv("data_ignored/secondary/winDf_SampleValues.csv")
sampleDf$indOg <- sampleDf$ind
sampleDf$ind <- gsub("phal_|-[[:digit:]]$|_1$","",sampleDf$indOg)
win_nonDivDf$id <- gsub("_0+","_",gsub("phal_|-[[:digit:]]$|_1$","",win_nonDivDf$Group.1))#gsub("_0+","_",gsub("phal_|-[[:digit:]]$|_1$","",win_nonDivDf$Group.1))
idCheck <- sort(win_nonDivDf$id[!(win_nonDivDf$id)%in%(sampleDf$ind)])
if(length(idCheck)!=0){
  stop("Id mismatch")
}else{
  win_nonDivDf$coverage <- sampleDf$newSamplePeakMean[match(win_nonDivDf$id,sampleDf$ind)]
  win_nonDivDf$madLogic <- win_nonDivDf$x>0.102
  plot(win_nonDivDf$coverage,win_nonDivDf$x,ylab="mad at non-divergent loci",xlab="peak coverage")
  text(win_nonDivDf$coverage[win_nonDivDf$madLogic],win_nonDivDf$x[win_nonDivDf$madLogic],
       labels = win_nonDivDf$id[win_nonDivDf$madLogic],col="red")
  points(win_nonDivDf$coverage,win_nonDivDf$x)
  print(tail(win_nonDivDf[order(win_nonDivDf$x),],20))
}


#### Making summary values ####
winDf_presub<- winDf
winDf_presub<- winDf_presub[order(winDf_presub$id,winDf_presub$chr,winDf_presub$min),]
#winDf_presub$div <- winDf_presub$median<min(cutoffs)|winDf_presub$median>max(cutoffs)

winDf_presub$divMinus2 <- c(F,F,winDf_presub$div[1:(nrow(winDf_presub)-2)])
winDf_presub$divMinus1 <- c(F,winDf_presub$div[1:(nrow(winDf_presub)-1)])
winDf_presub$divPlus1 <- c(winDf_presub$div[2:(nrow(winDf_presub))],F)
winDf_presub$divPlus2 <- c(winDf_presub$div[3:(nrow(winDf_presub))],F,F)

winDf_presub$Win3Sum <- winDf_presub$divMinus1 + winDf_presub$div +  winDf_presub$divPlus1
winDf_presub$Win5Sum <- winDf_presub$divMinus1 + winDf_presub$div +  winDf_presub$divPlus1 + winDf_presub$divMinus2 + winDf_presub$divPlus2

win3Max <- winregate(winDf_presub$Win3Sum,by=list(winDf_presub$id_seq),max)
out <- c(paste0("Of Id+Seq combos: ",sum(win3Max$x>=3)," (",round(mean(win3Max$x>=3)*100,2),"% of ",nrow(win3Max)," combos)"))
winDf_noScaffold <- winDf_presub[winDf_presub$chrLogic,]
win3Max <- winregate(winDf_noScaffold$Win3Sum,by=list(winDf_noScaffold$id_seq),max)
out <- c(out,paste0("Of Id+chr combos (no scaffolds): ",sum(win3Max$x>=3)," (",round(mean(win3Max$x>=3)*100,2),"% of ",nrow(win3Max)," combos)"))
win3Max <- winregate(winDf_noScaffold$Win3Sum,by=list(winDf_noScaffold$id),max)
out <- c(out,paste0("Of Ids: ",sum(win3Max$x>=3)," (",round(mean(win3Max$x>=3)*100,2),"% of ",nrow(win3Max)," combos)"))

win5Max <- winregate(winDf_presub$Win5Sum,by=list(winDf_presub$id_seq),max)
out <- c(out,paste0("Of Id+Seq combos: ",sum(win5Max$x>=5)," (",round(mean(win5Max$x>=5)*100,2),"% of ",nrow(win5Max)," combos)"))
winDf_noScaffold <- winDf_presub[winDf_presub$chrLogic,]
win5Max <- winregate(winDf_noScaffold$Win5Sum,by=list(winDf_noScaffold$id_seq),max)
out <- c(out,paste0("Of Id+chr combos (no scaffolds): ",sum(win5Max$x>=5)," (",round(mean(win5Max$x>=5)*100,2),"% of ",nrow(win5Max)," combos)"))
win5Max <- winregate(winDf_noScaffold$Win5Sum,by=list(winDf_noScaffold$id),max)
out <- c(out,paste0("Of Ids: ",sum(win5Max$x>=5)," (",round(mean(win5Max$x>=5)*100,2),"% of ",nrow(win5Max)," combos)"))

write(out,file = "data/plotsummaryvalues.csv")


#### Summarize values per chr and line
winDf$closestMed <- round(winDf$median)
winDf$chr_cnv <- paste0(format(winDf$chr_orig,justify = "right"),"_",format(winDf$closestMed))
winDf$ind__chr_cnv <- paste0(winDf$id,"__",winDf$chr_cnv)

cnvDf <- data.frame(table(winDf$ind__chr_cnv))
cnvDf$id     <- trimws(gsub("__.*","",cnvDf$Var1))
cnvDf$chr_cn <- trimws(gsub(".*__","",cnvDf$Var1))
cnvDf$chr    <- trimws(gsub("_.*","",cnvDf$chr_cn))
cnvDf$cn     <- as.numeric(gsub(".*_","",cnvDf$chr_cn))

uniChr_cn <- unique(cnvDf$chr_cn)
for(i in 1:length(uniChr_cn)){
  if(i==1){outDf<-NULL}
  currChr_cn <- uniChr_cn[i]
  currDf <- cnvDf[cnvDf$chr_cn==currChr_cn,]
  outDf<-rbind(outDf,currDf[which.max(currDf$Freq),])
}
outDf <- outDf[outDf$Freq>5,]
outDf$chrLogic <- grepl("Chr",outDf$chr)
outDf <- outDf[order(-outDf$chrLogic,outDf$chr,outDf$cn),]
outDf$url <- winDf$singleLineMultiChr[match(outDf$id,winDf$id)]
write.csv(outDf,"data/linesWithMostCNVs.csv")




