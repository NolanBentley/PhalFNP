#LOad packages
library(dplyr)
library(ggplot2)
library(treemapify)
library(reshape2)


#Setup
wd <- "~/Experiments/PhalFNP/"; setwd(wd)
coverageFile <- "data/sampleDf.csv"
variantFile  <- "https://utexas.box.com/shared/static/yht8ojfafsab3btuq6xe2chw8by9onfo.csv"

#Load data
samDf <- read.csv(coverageFile)
df1   <- read.csv(variantFile )

#Mutation burden
mean(aggregate(df1$id,by = list(df1$sample),length)$x)

#Store modified values
df1$id_orig     <-df1$id
df1$sample_orig <-df1$sample

##### Clean data ######
###### Add / modify values #####
#Clean sample id
df1$sample <- gsub("Fil","FIL",gsub("M2_","M2-",gsub("\\.","_",df1$sample)))
df1$sample <- gsub("-.$","",df1$sample)
df1$mother <- gsub("(FIL.._[[:digit:]]+)_.*","\\1",df1$sample)
df1$dosage <- gsub("(FIL..).*","\\1",df1$mother)
df1$motherNum <- as.numeric(gsub("FIL.._","",df1$mother))
df1$familyLetter<- gsub("_M2|^M2$","",mapply(FUN = gsub,pattern=paste0(df1$mother,"_"),replacement="",x=df1$sample))
df1$sortableSample <- gsub(" ",0,paste0(df1$dosage,"_",format(df1$motherNum),"_",format(df1$familyLetter,justify = "right")))


#Remove bad samples
# Odd ids: "FIL30_774_M2b"   "FIL30_347_M2dup" "FIL30_717_M2b" 
removedSamples <- c(
  "FIL30_390_M2", #high sd of depth, not in variants data
  "FIL30_400_M2", #high sd of depth, not in variants data
  "FIL30_591_M2", #high sd of depth, not in variants data
  "FIL30_716_M2", #high sd of depth, not in variants data
  "FIL30_347_M2dup", #redundent
  "FIL30_717_M2b", #redundent
  "FIL30_774_M2b" #redundent
)
unique(df1$sample[which(df1$sample%in%removedSamples)])


df1[df1$motherNum%in%c(774,347,717),]

#Add in analysis data 
df1$analysis <- "GATK"
df1$analysis[is.na(df1$QD)]<-"DELLY"

#Add in positions
df1$begPos  <- as.numeric(gsub("->.*","",df1$pos))
df1$endPos  <- as.numeric(gsub(".*->","",df1$pos))
df1$endPosAlt <- 0
df1$lengthAlt<-df1$length
df1$lengthAlt[is.na(df1$length)] <- (nchar(df1$ref_allele)-1)[is.na(df1$length)]
df1$endPosAlt[df1$type%in%c("DEL","DUP","INV","SNP","INDEL")]<-(df1$begPos+df1$lengthAlt)[df1$type%in%c("DEL","DUP","INV","SNP","INDEL")]
df1$endPosAlt[df1$type%in%c("INS")]<-(df1$begPos+df1$length)[df1$type%in%c("INS")]-1
df1$endPosAlt[df1$type%in%c("BND","INDEL,OTHER")]<-NA
df1$EndMinusBeg <- df1$endPos - df1$begPos
df1$posDiff  <- df1$endPosAlt-df1$endPos
df1$posMatch <- df1$posDiff==0

head(df1[which(df1$motherNum%in%c(347)&df1$analysis=="GATK"),c("chrom","pos","type","begPos","endPos","length","lengthAlt")],40)

#Clean ID value
df1$posId <- paste0(df1$chrom,"-",gsub(" ","0",paste0(format(df1$begPos),"_",format(df1$endPosAlt))))
df1$id    <- paste0(df1$posId,".",df1$alt_allele_id,"|",df1$id_orig)
df1$sample_locusId <- paste0(df1$sample,"_",df1$id)

#Reorder
df1 <- df1[order(df1$chrom,df1$begPos,df1$endPos,df1$sortableSample,df1$id),]

#Create column of non-redundant sample-id combos
duplicatedSampleLocusIds <- df1$sample_locusId[duplicated(df1$sample_locusId)]
df1$hasDupeSamLocId  <- df1$sample_locusId%in%duplicatedSampleLocusIds 
df1$nr <- !duplicated(df1$sample_locusId)

#Detect reoccurring variants
hasDupes <- function(x){x%in%(x[duplicated(x)])}
df1$id_duplicated<-hasDupes(df1$id)

#Detect near gene
df1$nearGene <- !(is.na(df1$gene)|df1$gene=="")|!is.na(df1$overlapping.genes)

#Test for chromosome localized
df1$chromLogic <- grepl("^Chr..$|Chr.*->Chr",x = df1$chrom)

#Make subset of just unique variant id + genotype combinations excluding scaffolds
df2<-df1[df1$nr&df1$chromLogic,]
length(unique(df2$id)) #Unique variant IDs
nrow(df2) #Unique sample ID + variant ID combos
df3 <- df2[df2$nearGene,]
length(unique(df3$id)) #Unique variant IDs near genes

#Load gene data
if(!file.exists("data/geneDf.csv")){source("scripts/functions/generateGeneDf.R")}
geneDf <- read.csv("data/geneDf.csv")

#Test gene proximity
df3$geneStart<-geneDf$start[match(df3$gene,geneDf$id)]
df3$geneEnd  <-geneDf$end[match(df3$gene,geneDf$id)]
plot((df3$begPos - df3$geneStart)[df3$analysis=="DELLY"],
     (df3$endPos - df3$geneEnd)[df3$analysis=="DELLY"],
     xlim=c(-3000,3000),ylim=c(-3000,3000))
abline(v = -1500)
abline(h = 1500)

#Make table GATK
gatkDf <- df3[df3$analysis=="GATK"&!duplicated(df3$id),]
gatkDf$tabledVals <- paste0(gatkDf$impact,"|",gatkDf$region)
gatkDf$tabledVals <- gsub("\\+.*","+OTHER",gatkDf$tabledVals)
gatkTable <- as.matrix(table(gatkDf$tabledVals,gatkDf$type))
gatkTable <- gsub(" ","",format(gatkTable,big.mark = ","))
gatkTableDf <- data.frame(
  Impact=gsub("\\|.*","",rownames(gatkTable)),
  Region=gsub(".*\\|","",rownames(gatkTable)),
  INDEL=as.vector(gatkTable[,"INDEL"]),
  SNP=as.vector(gatkTable[,"SNP"])
)
write.csv(gatkTableDf,file = "data/tableGatk_GenicVariantSummary.csv")

#Make tableDelly
dellyDf <- df3[df3$analysis=="DELLY"&!duplicated(df3$id),]
nrow(dellyDf) #Total unique variants in Delly analysis
sum(dellyDf$type=="BND") #Total number of transloaction
dellyDf$lengthBinNum <- 10^(ceiling(log10(dellyDf$length)))
dellyDf$lengthBin <- gsub(" ","",format(ceiling(dellyDf$lengthBinNum),big.mark = ",",scientific = F))
dellyDf$lengthBin <- paste0("<=",gsub(",000$","Kb",gsub(",000,000$","Mb",dellyDf$lengthBin)))
tableDelly <- gsub(" ","",format(as.matrix(table(dellyDf$type,dellyDf$lengthBin)),big.mark = ","))
tableDelly <- tableDelly[,unique(dellyDf$lengthBin[order(dellyDf$lengthBinNum)])]
write.csv(tableDelly,file = "data/tableDelly_GenicVariantSummary.csv")


####Make dosage plots####
#Calculate mutation count per sample
df2SamAgg2 <- aggregate(df2$sample,by=list(sample=df2$sample),length)
df2SamAgg2$dosage <- gsub("_.*","",df2SamAgg2$sample)
table(df2SamAgg2$dosage)#How many of each dosage
mean(df2SamAgg2$x[df2SamAgg2$dosage=="FIL20"])#Mean at dosage 20
sd(df2SamAgg2$x[df2SamAgg2$dosage=="FIL20"])
mean(df2SamAgg2$x[df2SamAgg2$dosage=="FIL30"])#Mean at dosage 30
sd(df2SamAgg2$x[df2SamAgg2$dosage=="FIL30"])

#Calculate mutation count per sample & analysis
df2SamAgg <- aggregate(df2$sample,by=list(sample=df2$sample,analysis=df2$analysis),length)
df2SamAgg$dosage <- gsub("_.*","",df2SamAgg$sample)

#Violin plot
ggplot(df2SamAgg,aes(analysis,x,fill=dosage))+
  geom_violin(draw_quantiles = c(0.5),alpha=0.5)+
  theme_bw()+
  labs(x="Analysis",y="Variants per genotype")
#Density
ggplot(df2[which(df2$lengthAlt>0),],aes(lengthAlt,fill=dosage))+
  geom_density(alpha=0.5)+
  theme_bw()+
  scale_x_continuous(trans="log10",
                     breaks = 10^(0:6),
                     minor_breaks = rep(0:5*2,5)*10^sort(rep(0:4,6)))+
  labs(x="Variant length (bp)",y="Density of variants longer than 1 bp per dosage")





#Tree plot of variant types
treeFun<-function(x){
  colnames(x)[colnames(x)=="value"]<-"count"
  p1 <- ggplot(x, aes(area = count, fill = count, label=label)) +
    geom_treemap()+
    geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                      grow = F)+
    scale_fill_viridis_c(trans="log10",direction = -1)
  return(p1)
}
treewidth <- 8
treeratio <- 3/5 
tab3 <- melt(table(type=df2$type,impact=df2$impact))
plottedDf <- tab3[tab3$value>0,]
plottedDf$label <- paste0("Type: ",plottedDf$type,"\n","Impact: ",plottedDf$impact,"\nn=",gsub(" ","",format(plottedDf$value,big.mark = ",")))
ggsave(treeFun(x = plottedDf),filename = "data/treemap.png",
       width = treewidth,height = treewidth*treeratio,
       dpi = 600,units = "in")

#
uniChr <- unique(df1$chrom[grepl("^.hr..$",df1$chrom)])
for(i in uniChr){
  for(j in uniChr){
    bndMat[i,j]<-sum(df1$chrom==paste0(i,"->",j))
  }
}