##### Setup environent ####
#wd
wd <- "~/Experiments/PhalFNP/"
coverageFile <- "data_ignored/secondary/aggDf_SampleValues.csv"
setwd(wd)

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
mainText <- paste0("\n",
  sum(samDf$cp_logic)," (",round(mean(samDf$cp_logic)*100,3),"%) of ",nrow(samDf),
  " genotypes are likely cross-pollinated"
)

#Establish the parent of likely self-s
hetPropCutoff2 <- 0.9
samDf$father <- NA
samDf$likelySelf <- samDf$lfm_n>=(totCntCutoff*2)&samDf$lfm_hetProp<=hetPropCutoff2
samDf$father[samDf$likelySelf] <- samDf$mother[samDf$likelySelf]

#Add in coverage analysis
samCoverage <- read.csv(coverageFile)
samCoverage$indMod <- gsub("-","_",gsub("phal_","",samCoverage$ind))
samDf$originMod <- gsub("-","_",gsub("Fil","FIL",gsub("\\.","-",samDf$origin)))
head(sort(samCoverage$indMod[!samCoverage$indMod%in%samDf$originMod]))
head(sort(samDf$originMod[!samDf$originMod%in%samCoverage$indMod]))
samCoverage$match <- match(samCoverage$indMod,samDf$originMod)
samDf$coverage <- NA
samDf$coverage[samCoverage$match[!is.na(samCoverage$match)]] <- samCoverage$newSamplePeakMean[!is.na(samCoverage$match)]

#Make a graphic
library(ggExtra)
library(ggplot2)

samDf$cp_status <- "Ambiguous"
samDf$cp_status[samDf$likelySelf]<-"Self"
samDf$cp_status[samDf$cp_logic  ]<-"Cross"
p1 <- ggplot(samDf[order(-samDf$coverage),],
             aes(lfm_hetProp,lfm_n,fill=cp_status,size=coverage))+
  #theme_bw()+#
  labs(x="Proportion of filtered variants that are heterozygous",
       y=paste0("Number of filtered variants"),
       fill="Parentage estimate:",
       size="Peak coverage:",
       title=mainText)+
  geom_vline(xintercept = c(hetPropCutoff2,hetPropCutoff),color=c("forestgreen","purple4"))+
  geom_hline(yintercept = c(totCntCutoff*2,totCntCutoff),color=c("forestgreen","purple4"))+
  geom_point(color="white",fill="white",mapping = aes(size=coverage+0.2*max(coverage,na.rm=T)))+
  geom_point(alpha=0.3,shape=21)+
  scale_size(range = c(0.5,8))+
  scale_y_continuous(breaks = seq(0,200,by=20))+
  scale_x_continuous(breaks = seq(0,1,by=0.1))+
  theme(legend.position="bottom",
        legend.box="vertical", 
        legend.margin=margin(),
        plot.title=element_text(hjust = 0.5,size = 10))
p2 <- ggMarginal(p1,type="histogram",size=5,bins=100);p1
ggsave("data/heterozygosity.png",p2,width = 88*2,height = 88*2.2,units = "mm",dpi = 600)


#Calculate similarity values across nr / lfm variants
relaDf <- similarityFun(df1[df1$nr&df1$lowFreqMutation,])
relaDf$logicSum <-rowSums(relaDf[,grep("^(het|homo)_(het|homo)$",colnames(relaDf))])
relaDf$totSum   <-rowSums(relaDf[,grep("tot_(i|j)",colnames(relaDf))])
relaDf$totMean  <- relaDf$totSum/2
relaDf$totProp  <- relaDf$logicSum/relaDf$totMean

#Make smaller
hist(relaDf$totProp,ylim=c(0,1000),10000);abline(v = 0.05,col="red")
relaDf <- relaDf[relaDf$totProp>0.05,]

#Filter and organize
relaDf <- relaDf[relaDf$iSam!=relaDf$jSam,]
relaDf <- relaDf[order(-relaDf$totProp),]

#Cross to inbred analysis
chiFromTot <- function(rel,totCol,divVector){
  exp_hethet<-totCol/(divVector[1])
  exp_hethom<-totCol/(divVector[2])
  exp_homhet<-totCol/(divVector[3])
  exp_homhom<-totCol/(divVector[4])
  chiSum <- rowSums(na.rm = T,cbind(
    (rel$het_het  - exp_hethet)^2/exp_hethet,
    (rel$het_homo - exp_hethom)^2/exp_hethom,
    (rel$homo_het - exp_homhet)^2/exp_homhet,
    (rel$homo_homo- exp_homhom)^2/exp_homhom
  ))
  return(chiSum)
}
table(samDf$cp_status)
  
relaDf$crToIb <- chiFromTot(relaDf,relaDf$tot_j,c(3,6,0,0))
relaDf$ibToCr <- chiFromTot(relaDf,relaDf$tot_i,c(3,0,6,0))
relaDf$ibToIb <- chiFromTot(relaDf,(relaDf$tot_j+relaDf$tot_i)/2,c(3,6,6,12))
relaDf$fs     <- chiFromTot(relaDf,relaDf$tot_i,c(2,0,0,0))+
  chiFromTot(relaDf,relaDf$tot_j,c(2,0,0,0))
relaDf$minTest <- apply(cbind(relaDf$crToIb,relaDf$ibToCr,relaDf$ibToIb,relaDf$fs),1,min,na.rm=T)

#Add in known parentage
relaDf$iSam_mother <- samDf$mother[match(relaDf$iSam,samDf$sample)]
relaDf$iSam_father <- samDf$father[match(relaDf$iSam,samDf$sample)]
relaDf$jSam_mother <- samDf$mother[match(relaDf$jSam,samDf$sample)]
relaDf$jSam_father <- samDf$father[match(relaDf$jSam,samDf$sample)]
relaDf <- data.frame(MatSib=relaDf$iSam_mother==relaDf$jSam_mother,relaDf[,colnames(relaDf)!="MatSib"])

manualAssignment <- cbind(relaDf,"...",samDf[match(relaDf$iSam,samDf$sample),],"...",samDf[match(relaDf$jSam,samDf$sample),])
write.csv(manualAssignment,file = "data/highRelas.csv")
manualAssignment <- manualAssignment[is.na(manualAssignment$iSam_father),]
write.csv(manualAssignment,file = "data/fatherlessRelas.csv")

#Taking notes on relationships
write.csv(samDf,file = "data/sampleDf.csv")
