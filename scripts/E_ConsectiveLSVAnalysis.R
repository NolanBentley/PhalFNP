#Setup environment
wd <- "~/Experiments/PhalFNP";setwd(wd)
kVal <- 7
divCutoffs <- c(1.5,2.5)
outDir  <- "./depthImages"
unzippedCovPath <- "./data_ignored/secondary/cnv/unzipped"
winName <- "7x10kb"
badLocusSpacing <- 2500

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

#Load data
compDV2 <- read.csv("data_ignored/secondary/compDV2.csv")

#Check for k consecutive divergent intervals
xMax <- ceiling(max(compDV2$mid)/1000000)
binTable <- table(compDV2$binExt)
dir.create("data/consecutiveLSVPlots",showWarnings = F)
for(kVal in seq(3,25,by=2)){
  #
  compDV2$kIntLogic <- compDV2$binExt%in%as.numeric(names(which(binTable>=kVal)))
  #Full "divergent" plot
  set.seed(1)
  plottedDf <- compDV2[compDV2$kIntLogic,]
  plottedDf <- plottedDf[sample(1:nrow(plottedDf)),]
  p1 <- ggplot(plottedDf,aes(mid/1000000,CN,color=id))+
    geom_hline(yintercept = 2,color="red")+
    geom_point()+
    theme_bw()+
    facet_wrap(~chr,ncol=3)+
    scale_y_continuous(limits = c(0,10  ),breaks = seq(0,10))+
    scale_x_continuous(limits = c(0,xMax),breaks = seq(0,xMax,by=5))+
    scale_color_discrete(guide="none")+
    labs(y=paste0("CN (minimum of ",kVal, " intervals = ",format(nrow(plottedDf),big.mark = ","),")"),
         x="Position in Mb")+
    geom_rug(data = geneDf,mapping = aes(x=(end+start)/2,y=NULL,color=NULL),inherit.aes = F)
  ggsave(p1,file=paste0("data/consecutiveLSVPlots/fullCnvPlot_",kVal,".png"),width = 10,height = 8,dpi = 600)
  print(nrow(plottedDf))
}
