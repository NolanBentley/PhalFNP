##### Setup environent ####
#wd
wd <- "~/Experiments/PhalFNP/"
coverageFile <- "data_ignored/secondary/aggDf_SampleValues.csv"
setwd(wd)

#Load packages
library(ggExtra)
library(ggplot2)
library(plotly)
library(treemap)

#Load data
source("./scripts/functions/similarityCalc.R")
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

#Clean variant id
df1$analysis <- "GATK"
df1$analysis[is.na(df1$QD)]<-"DELLY"
df1$id      <- paste0(df1$chrom,"-",df1$pos,"_",df1$ref_allele,"->",df1$alt_allele)
df1$sample_locusId <- paste0(df1$sample,"_",df1$id)

#Add in positions
df1$begPos  <- as.numeric(gsub("->.*","",df1$pos))
df1$endPos  <- as.numeric(gsub(".*->","",df1$pos))

#Reorder
df1 <- df1[order(df1$chrom,df1$begPos,df1$endPos,df1$sample,df1$id),]

#Only keep one row describing each sample_mutation-locus combination
duplicatedSampleLocusIds <- df1$sample_locusId[duplicated(df1$sample_locusId)]
df1$hasDuplicateLocusId  <- df1$sample_locusId%in%duplicatedSampleLocusIds 
df1$nr <- !duplicated(df1$sample_locusId)

#Detect reoccurring variants
df1$id_duplicated<-hasDupes(df1$id)

