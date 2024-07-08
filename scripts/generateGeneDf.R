#Setup
wd <- "~/Experiments/PhalFNP/"
gff3File <- "data_ignored/primary/assembly/Phallii_590_v3.2.gene.gff3"
attributesRegex <- "^ID\\=(.*?)\\.v3.2\\;.*$"
phytozomePrefix <- "https://phytozome-next.jgi.doe.gov/report/transcript/Phallii_v3_2/"
outFile <- "data/geneDf.csv"

#Make file
setwd(wd)
geneDf <- read.delim(gff3File,header = F,skip = 3)
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
geneDf$id <- gsub(attributesRegex,"\\1",geneDf$attributes)
geneDf$link <- paste0(phytozomePrefix,geneDf$id,".1")
write.csv(geneDf,outFile)