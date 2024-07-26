#Set variables
wd <- "~/Experiments/PhalFNP/"
minInitFilter<-0.001
maxInitFilter<-0.02
vcfFile    <- "data_ignored/primary/vcf/genotype_gvcfs.vcf"

#Setup
setwd(wd)
library(data.table)

#Figure out heading
header<-readLines(vcfFile,n = 2000)
firstNonComment<-which(!grepl("^#",header))[1]#header[330:332]
nCol <-length(strsplit(header[firstNonComment-1],split = "\t")[[1]])
nSam <- nCol-9

#Initial filter
minInit <- ceiling(minInitFilter*nSam)
maxInit <- floor(maxInitFilter*nSam)
vcfSubFile <- paste0(gsub("\\.vcf","",vcfFile),"_gte",minInit,"xlte",maxInit,"mac.vcf")
initFilterCmd <- paste0(
  "bcftools view --min-ac ",minInit,
  " --max-ac ",maxInit,
  " ",vcfFile," -o ",vcfSubFile
)
if(!file.exists(vcfSubFile)){system(initFilterCmd,wait = T)}

#Load data
vcf <- fread(vcfSubFile,sep = "\t",header = T,skip = firstNonComment-1)
colnames(vcf)

#Clean column names
colnames(vcf)<- gsub("^#|_1$|-[1-4]$","",gsub("Fil","FIL",gsub("\\.","_",colnames(vcf))))

#Isolate and clean genos
colVals <- which(colnames(vcf)=="FIL2"):ncol(vcf)
genos   <- apply(vcf[,..colVals],2,gsub,pattern="\\:.*",replacement="")

#Calculate locus variables
hetLogic  <- genos == "0/1"
homoLogic <- genos == "1/1"
missLogic <- genos == "./."
minorAlleleCnt <- (hetLogic+homoLogic*2)/(!missLogic)
vcf$hetFreq    <- rowMeans(hetLogic )
vcf$homoMAFreq <- rowMeans(homoLogic)
vcf$missFreq   <- rowMeans(missLogic)
vcf$maf        <- rowMeans(minorAlleleCnt,na.rm=T)
hist(vcf$maf,1000)
hist(vcf$missFreq,ylim=c(0,1000),1000)

#Isolate subset of mid maf alleles
vcf$begPos <- vcf$POS-nchar(vcf$REF)+1
vcf$endPos <- vcf$POS+nchar(vcf$REF)-1
vcf$marker <-paste0(vcf$CHROM,"_",vcf$begPos,"-",vcf$endPos)
if(any(duplicated(vcf$marker))){
  stop("Need to make marker names unique")
}else{
  rownames(vcf)<-vcf$marker
}

head(vcf$marker)
head(vcf$POS)
head(vcf$REF)
head(vcf$ALT)
vcf$sub1 <- 
  vcf$maf>=minInitFilter&
  vcf$maf<=maxInitFilter&
  genos[,"FIL2"]=="0/0"&
  vcf$missFreq<0.1
vcf2   <- as.data.frame(vcf[vcf$sub1,])
genos2 <- genos[vcf$sub1,]
mac2   <- minorAlleleCnt[vcf$sub1,]
rownames(vcf2)  <- rownames(vcf)[vcf$sub1]
rownames(genos2)<-rownames(vcf2)
rownames(mac2  )<-rownames(vcf2)

#write files
saveRDS(vcf2  ,"./data_ignored/secondary/midFreqVariants.rds")
saveRDS(genos2,"./data_ignored/secondary/midFreqGenos.rds")
saveRDS(mac2  ,"./data_ignored/secondary/midFreqMACs.rds")

