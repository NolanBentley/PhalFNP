##### Setup environent ####
#wd
wd <- "~/Experiments/PhalFNP/"
setwd(wd)

#Load data
source("./scripts/similarityCalc.R")
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

#For checking
df1$codedSample <- paste(df1$mother,df1$sample,df1$sample_orig,sep = "|")
samDf <- data.frame(mother=df1$mother,sample=df1$sample,origin=df1$sample_orig,code=df1$codedSample)
samDf <- samDf[!duplicated(samDf$code),]
samDf <- samDf[order(samDf$mother,samDf$sample,samDf$orig),]
samDf$n_tot      <- cntOfYInX(df1$codedSample,samDf$code)
samDf$n_delly    <- cntOfYInX(df1$codedSample[df1$analysis=="DELLY"],samDf$code)  
samDf$n_gatk     <- cntOfYInX(df1$codedSample[df1$analysis=="GATK"],samDf$code)  
samDf$nr_n       <- cntOfYInX(df1$codedSample[df1$nr                     ],samDf$code)
samDf$nr_het     <- cntOfYInX(df1$codedSample[df1$nr & df1$var_count==2  ],samDf$code)
samDf$nr_homo    <- cntOfYInX(df1$codedSample[df1$nr & df1$var_count==1  ],samDf$code)  
samDf$nr_uni     <- cntOfYInX(df1$codedSample[df1$nr & !df1$id_duplicated],samDf$code) 
samDf$nr_uniHet  <- cntOfYInX(df1$codedSample[df1$nr & !df1$id_duplicated & df1$var_count==1],samDf$code)  
samDf$nr_uniHomo <- cntOfYInX(df1$codedSample[df1$nr & !df1$id_duplicated & df1$var_count==2],samDf$code)
samDf$motherDL <- hasDupes(samDf$mother)
samDf$sampleDL <- hasDupes(samDf$sample)

#Find high-freq mutations
id_table <-table(df1$id[df1$nr])
hist(id_table,10000,ylim=c(0, 10),xlim=c(0,100))
abline(v = 9.5,col="red")
df1$lowFreqMutation <- df1$id%in%names(id_table)[id_table<9.5]

#Find likely outcrosses based on over representation of heterozygous genotypes
dfCP <- df1[df1$lowFreqMutation & df1$nr & df1$type%in%c("SNP","INDEL"),]

samDf$lfm_n       <- cntOfYInX(dfCP$codedSample,samDf$code)
samDf$lfm_het     <- cntOfYInX(dfCP$codedSample[dfCP$var_count==1],samDf$code)
samDf$lfm_homo    <- cntOfYInX(dfCP$codedSample[dfCP$var_count==2],samDf$code)
samDf$lfm_hetProp <- samDf$lfm_het/samDf$lfm_n
hetPropCutoff <- 0.999
totCntCutoff  <- 10
samDf$cp_logic <- samDf$lfm_hetProp>=hetPropCutoff & samDf$lfm_n>=totCntCutoff
df1$cp_logic   <- df1$sample%in%(samDf$sample[which(samDf$cp_logic)])
table(samDf$cp_logic,useNA = "ifany")
mainText <- paste0(
  "Of the samples, ",
  sum(samDf$cp_logic)," (",round(mean(samDf$cp_logic)*100,3),"%) of ",nrow(samDf),
  " genotypes are likely cross-pollinated"
)

#Make a graphic
library(ggExtra)
library(ggplot2)
p1 <- ggplot(samDf,aes(lfm_hetProp,lfm_n,color=cp_logic))+
  theme_bw()+#
  labs(x="Proportion of filtered variants heterozygous",
       y=paste0("Number of variants after filtering to low frequency mutations (lfm)"),
       color="Likely cross-pollinated:",
       title=mainText)+
  geom_vline(xintercept = hetPropCutoff,color="forestgreen")+
  geom_hline(yintercept = totCntCutoff,color="purple4")+
  geom_point(color="white",size=3)+
  geom_point(alpha=0.5,size=2.5)+
  theme(legend.position = "bottom")+
  scale_y_continuous(breaks = seq(0,200,by=20))+
  scale_x_continuous(breaks = seq(0,1,by=0.1))
p2 <- ggMarginal(p1,type="histogram",size=5,bins=100)
ggsave("data/heterozygosity.png",p2,width = 88*2,height = 88*2.2,units = "mm",dpi = 600)


#Establish the parent of likely self-s
hetPropCutoff2 <- 0.9
sampleAgg$mother <- df1$mother[match(sampleAgg$Group.1,df1$sample)]
sampleAgg$father <- NA
sampleAgg$likelySelf <- sampleAgg$totCnt>=(totCntCutoff*2)&sampleAgg$hetProp<=hetPropCutoff2
sampleAgg$father[sampleAgg$likelySelf] <- sampleAgg$mother[sampleAgg$likelySelf]

plot(sampleAgg$hetProp,sampleAgg$homoCnt*2-sampleAgg$hetProp,col=c("black","red")[sampleAgg$likelySelf+1])

#Manually set known relationships based on best-judgement
table(sampleAgg$father,useNA = "ifany")

uniFatherless <- sampleAgg$Group.1[is.na(sampleAgg$father)]
df1Homo <- df1[df1$var_count==2,]
for(i in 1:length(uniFatherless)){
  if(i==1){
    basDf <- NULL
  }
  currVars  <- df1[df1$sample==uniFatherless[i],]
  currAgg     <- aggregate(df1$id    ,by=list(df1$sample    ),function(x,y){sum(x%in%y)},y=currVars$id)
  currHomoAgg <- aggregate(df1Homo$id,by=list(df1Homo$sample),function(x,y){sum(x%in%y)},y=currVars$id)
  currAgg$homo_x <- 0
  currAgg$homo_x [match(currHomoAgg$Group.1,currAgg$Group.1)]<- currHomoAgg$x
  currAgg$fatherlessName   <- currVars$sample[1]
  basDf<- rbind(basDf,currAgg)
  print(i)
}
basDf <- basDf[basDf$Group.1!=basDf$fatherlessName,]
basDf <- basDf[basDf$x>0,]

basDf$flMother <- df1$mother[match(basDf$fatherlessName,df1$sample)]
basDf$pfMother <- df1$mother[match(basDf$Group.1       ,df1$sample)]
basDf$pfFather <- sampleAgg$father[match(basDf$Group.1       ,sampleAgg$Group.1)]
basDf$n12      <- sampleAgg$totCnt[match(basDf$fatherlessName,sampleAgg$Group.1)]
basDf$n11      <- sampleAgg$totCnt[match(basDf$Group.1       ,sampleAgg$Group.1)]

basDf$inbredToSib <- (basDf$n11*0.5-basDf$x)^2/(basDf$n11*0.5)
hist(basDf$inbredToSib,1000)


basDf$n1 <- basDf$n11/0.75
basDf$n2 <- (basDf$n12 - 0.5*(basDf$n1))/0.5

basDf$n2 - basDf$n1

# P1xP2  vs. P1xP1 assuming P1 fully heterozygous and P1 and P2 have similar frequencies of mutation
#P1xP2 heterozygous at half of P1 
#P1xP2 heterozygous at half of P1 

#P1xP1 homozygous mutant at 33% of P1
#P1xP1 heterozygous at 

View(basDf[order(-basDf$x),])

basDf$x

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
