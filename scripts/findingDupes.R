#wd
wd <- "~/Experiments/PhalFNP/"
setwd(wd)

#Load data
source("./scripts/similarityCalc.R")
variantFile <- "https://utexas.box.com/shared/static/yht8ojfafsab3btuq6xe2chw8by9onfo.csv"
df1 <- read.csv(variantFile)

#Save modified columns
if(is.null(df1$id_orig)){df1$id_orig<-df1$id}

#Clean up data
df1$analysis <- "GATK"
df1$analysis[is.na(df1$QD)]<-"DELLY"
df1$idFull  <- paste0(df1$id_orig,":",df1$chrom,"-",df1$pos,"_",df1$ref_allele,"->",df1$alt_allele,"_",df1$analysis,"_",df1$distance,"_",df1$Transcript_ID)
df1$id      <- paste0(df1$chrom,"-",df1$pos,"_",df1$ref_allele,"->",df1$alt_allele)
df1$sample_locusId <- paste0(df1$sample,"_",df1$id)
df1$sample_idFull  <- paste0(df1$sample,"_",df1$idFull)
df1$begPos  <- as.numeric(gsub("->.*","",df1$pos))
df1$endPos  <- as.numeric(gsub(".*->","",df1$pos))
df1 <- df1[order(df1$chrom,df1$begPos,df1$endPos,df1$sample,df1$id),]

#Only keep one row describing each sample_mutation-locus combination
df1<- df1[!duplicated(df1$sample_locusId),]

#Remove high freq mutations
id_table <-table(df1$id)
hist(id_table,10000,ylim=c(0, 100),xlim=c(0,100))
abline(v = 30.5,col="red")
df1$lowFreqMutation <- df1$id%in%names(id_table)[id_table<30.5]
df1<- df1[df1$lowFreqMutation,]

#Remove samples without any overlaps
df1$id_hasDupe <- df1$id%in%(df1$id[duplicated(df1$id)])
df1$sample_inDupe <- df1$sample%in%(df1$sample[df1$id_hasDupe])
mean(df1$sample_inDupe)
#df1 <- df1[df1$sample_inDupe,]


#Clean sample id
sort(unique(df1$sample))
df1$sample <- gsub("Fil","FIL",gsub("M2_","M2-",gsub("\\.","_",df1$sample)))
if(all(table(gsub("-.*","",df1$sample))==table(df1$sample))){
  df1$sample<-gsub("-.*","",df1$sample)
}
df1$mother <- gsub("(FIL.._[[:digit:]]+?)_.*","\\1",df1$sample)


#Find likely outcrosses based on over representation of heterozygous genotypes
df0 <- df1
df0 <- df0[df0$type%in%c("SNP","INDEL"),]
sampleAgg         <- aggregate(df0$var_count-1,by=list(df0$sample),mean)
sampleAgg$homoCnt <- aggregate(df0$var_count-1,by=list(df0$sample),sum )$x
sampleAgg$hetCnt  <- aggregate(2-df0$var_count,by=list(df0$sample),sum )$x
sampleAgg$totCnt  <-sampleAgg$homoCnt+sampleAgg$hetCnt
sampleAgg$hetProp <- (sampleAgg$hetCnt)/sampleAgg$totCnt
hetPropCutoff <- 0.999
totCntCutoff  <- 10
df1$outcross <- df1$sample%in%(sampleAgg$Group.1[sampleAgg$hetProp>=hetPropCutoff&sampleAgg$totCnt>=totCntCutoff])
hist(sampleAgg$hetProp,1000,main=paste0(
  sum(sampleAgg$hetProp>hetPropCutoff),
  " (",round(mean(sampleAgg$hetProp>hetPropCutoff)*100,3),"%) of ",nrow(sampleAgg),
  " genotypes are likely outcrossed")
)
abline(v = hetPropCutoff,col="pink",lty=2)
plot(sampleAgg$hetProp,sampleAgg$totCnt)
abline(v = hetPropCutoff,col="pink",lty=2)
abline(h = totCntCutoff,col="pink",lty=2)

#Establish the parent of likely self-s
hetPropCutoff2 <- 0.9
sampleAgg$mother <- df1$mother[match(sampleAgg$Group.1,df1$sample)]
sampleAgg$father <- NA
sampleAgg$likelySelf <- sampleAgg$totCnt>=(totCntCutoff*2)&sampleAgg$hetProp<=hetPropCutoff2
sampleAgg$father[sampleAgg$likelySelf] <- sampleAgg$mother[sampleAgg$likelySelf]

#Manually set known relationships based on best-judgement
sampleAgg$father

uniFatherless <- sampleAgg$Group.1[is.na(sampleAgg$father)]
for(i in 1:legnth(uniFatherless)){
  uniFatherless
}

#Look for 
df1VarTable <-aggregate(df1$var_count,by=list(df1$sample),table)
df1VarFreqs <- do.call("rbind",lapply(df1VarTable$x,function(x){c(het=x["1"],homo=x["2"])}))
df1VarFreqs[is.na(df1VarFreqs)]<-0
plot(df1VarFreqs[,1],df1VarFreqs[,2],col=c("black","red")[(df1VarTable$Group.1%in%(df1$sample[df1$outcross]))+1])

heatmap(samMat[samplesInvestigated,samplesInvestigated],scale = "none")

#Explore the results
diag(samMat)<-NA

highMax <-apply(samMat,1,max,na.rm=T)
highMaxMat <- samMat[highMax>0.6,highMax>0.6]
heatmap(highMaxMat,scale="none")

install.packages("heatmaply")
install.packages("htmlwidgets")
library(heatmaply)
library(htmlwidgets)
pMap <- heatmaply(highMaxMat)
pMap
