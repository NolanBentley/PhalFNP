##### Setup environment ####
#wd
wd <- "~/GitHub/PhalFNP/"
coverageFile <- "data_ignored/secondary/aggDf_SampleValues.csv"
setwd(wd)

#Load packages
library(ggplot2)
library(ggExtra)#install.packages("ggExtra")
library(plotly)
library(treemap)#install.packages("treemap")
library(ggvenn)#devtools::install_github("yanlinlin82/ggvenn")

#Load data
variantFile <- "https://utexas.box.com/shared/static/yht8ojfafsab3btuq6xe2chw8by9onfo.csv"
df1 <- read.csv(variantFile)

#Store modified values
df1$id_orig     <-df1$id
df1$sample_orig <-df1$sample

##### Clean data ######
###### Add / modify values #####
#Clean sample id
df1$sample <- gsub("Fil","FIL",gsub("M2_","M2-",gsub("\\.","_",df1$sample)))
df1$sample <- gsub("-.$","",df1$sample)
df1$mother <- gsub("(FIL.._[[:digit:]]+)_.*","\\1",df1$sample)

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

#Clean ID value
df1$posId <- paste0(df1$chrom,"-",gsub(" ","0",paste0(format(df1$begPos),"_",format(df1$endPosAlt))))
df1$id    <- paste0(df1$posId,".",df1$alt_allele_id,"|",df1$id_orig)
df1$sample_locusId <- paste0(df1$sample,"_",df1$id)

#Reorder
df1 <- df1[order(df1$chrom,df1$begPos,df1$endPos,df1$sample,df1$id),]

#Create column of non-redundant sample-id combos
duplicatedSampleLocusIds <- df1$sample_locusId[duplicated(df1$sample_locusId)]
df1$hasDupeSamLocId  <- df1$sample_locusId%in%duplicatedSampleLocusIds 
df1$nr <- !duplicated(df1$sample_locusId)

#Detect reoccurring variants
hasDupes <- function(x){x%in%(x[duplicated(x)])}
df1$id_duplicated<-hasDupes(df1$id)

#Detect near gene
df1$nearGene <- !(is.na(df1$gene)|df1$gene=="")

#Create tibble for plotting
tib1 <- as_tibble(df1[df1$nr,])
tib1$varInGatk  <- tib1$id%in%(tib1$id[tib1$analysis=="GATK" ])
tib1$varInDelly <- tib1$id%in%(tib1$id[tib1$analysis=="DELLY"])
tib1$"Var near gene" <- tib1$gene!=""&!is.na(tib1$gene)
tib1$"Impact not NA, Low, or Modifier" <- !tib1$impact%in%c(NA,"MODIFIER", "LOW")

ggplot(tib1, aes(
    A = `Var near gene`,
    B = `Impact not NA, Low, or Modifier`)
    ) +
    geom_venn() + theme_void() + coord_fixed()

exp <- tib1[!tib1$varInGatk&!tib1$"Impact not NA, Low, or Modifier",]
as.data.frame(exp[exp$type=="DEL"&exp$length>900000,])

ggplot(tib1, aes(
    A = varInGatk, 
    B = hasDupeSamLocId,
    C = varInDelly)
) +
    geom_venn() + theme_void() + coord_fixed()

hist(df1$lengthAlt[df1$analysis=="GATK"&df1$type=="INDEL"],1000,)


#Explorations
View(df1[df1$lengthAlt>20&df1$analysis=="GATK",])

df1[df1$lengthAlt==50&df1$nearGene&df1$analysis=="GATK",]

splitOverlaps0 <- strsplit(df1$overlapping.genes[!is.na(df1$overlapping.genes)],split = ",")
splitGenes <- strsplit(df1$gene,split = ",|\\+")

#Cleanup 
splitOverlaps  <- lapply(splitOverlaps0,function(x){sort(unique(gsub(".v3.2","",gsub("\\(.*","",x))))})

#Get unique
uniOverlap <- sort(unique(unlist(splitOverlaps)))
uniGenes   <- sort(unique(unlist(splitGenes)))

#Check for consistency with annotation
dfGene <- read.csv("data/geneDf.csv")
if(length(uniOverlap[!uniOverlap%in%dfGene$id])>0){stop("Id error!")}
if(length(uniGenes  [!uniGenes  %in%dfGene$id])>0){stop("Id error!")}

#Make summaries
dfMessage0 <- data.frame(
    message = c("Total unique mutations and sample combinations: ",
                "Total unique loci and measured allele combinations: ",
                "Total unique loci: "
    ))
dfMessage0$mea1 <- format(c(
    sum(df1$nr),
    length(unique(df1$id)),
    length(unique(df1$posId))
),big.mark = ",")
dfMessage0$mea2 <- format(c(
    sum(df1$nr[df1$nearGene|!is.na(df1$overlapping.genes)]),
    length(unique(df1$id[df1$nearGene|!is.na(df1$overlapping.genes)])),
    length(unique(df1$posId[df1$nearGene|!is.na(df1$overlapping.genes)]))
),big.mark = ",")
cat(paste0(dfMessage0$message,dfMessage0$mea1," (",dfMessage0$mea2," near genes)",collapse = "\n"))


dfMessage1 <-data.frame(message=c(
    "Number of unique annotations in overlap column: ",
    "Number of unique annotations in genes column: ",
    "Number added by the overlap column: "))
dfMessage1$mea1 <- format(c(
    length(uniOverlap),
    length(uniGenes),
    length(uniOverlap)- sum(uniGenes%in%uniOverlap)
),big.mark = ",")
cat(paste0(dfMessage1$message,dfMessage1$mea1,collapse = "\n"))


df1$lengthOfOverlap <- 0
df1$lengthOfOverlap[!is.na(df1$overlapping.genes)] <- unlist(lapply(splitOverlaps,length))

p1 <- ggplot(df1[df1$lengthOfOverlap>0&df1$nr,],aes(lengthAlt/1000,lengthOfOverlap))+
    geom_bin2d(bins = 40)+scale_fill_viridis_c(trans="log10")+
    labs(x="Length of variant in Kb",y="Number of genes in region")+
    theme_bw()+
    facet_wrap(~type)
ggsave(p1,file = "data/variantLengthVsNumberOfGenes.png",
       width = 8,height = 8 / 1.77,dpi=600)

df1[df1$type=="INS",]

df1$Analysis_Type <- paste0(df1$analysis,"-",df1$type)
p2 <- ggplot(df1[df1$nr,],aes(lengthAlt/1000,fill=as.factor(var_count)))+
    geom_histogram(bins=40)+
    labs(x="Length of variant in Kb",
         y="Number of variant rows",
         fill="Zygosity:\n1= Het\n2 = Homo")+
    theme_bw()+
    scale_y_continuous(trans="log10")+
    facet_wrap(~Analysis_Type)
ggsave(p2,file = "data/lengthVsTypeVsZygosity.png",
       width = 8,height = 8 / 1.77,dpi=600)
