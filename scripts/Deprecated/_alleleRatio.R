#Set variables
wd <- "~/Experiments/PhalFNP/"
vcfFile    <- "data_ignored/primary/vcf/genotype_gvcfs_gte1xlte15mac.vcf"

#Setup
setwd(wd)
library(data.table)

#Figure out heading
header<-readLines(vcfFile,n = 2000)
firstNonComment<-which(!grepl("^#",header))[1]#header[330:332]
nCol <-length(strsplit(header[firstNonComment-1],split = "\t")[[1]])
nSam <- nCol-9

#Initial filter
vcf <- fread(vcfFile,sep = "\t",header = T,skip = firstNonComment-2)
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

vcf$marker <- paste0(vcf$`#CHROM`,"_",vcf$POS)

hetLs<-list()
hetLs[[1]] <- vcf$`FIL30_155_M2-3`[grep("^0/1",vcf$`FIL30_155_M2-3`)]
names(hetLs[[1]]) <- vcf$marker[grep("^0/1",vcf$`FIL30_155_M2-3`)]
hetLs[[2]] <- vcf$`FIL30_154_M2-3`[grep("^0/1",vcf$`FIL30_154_M2-3`)]
hetLs[[3]] <- vcf$`FIL30_153_M2-3`[grep("^0/1",vcf$`FIL30_153_M2-3`)]
hetLs[[4]] <- vcf$`FIL30_157_M2-3`[grep("^0/1",vcf$`FIL30_157_M2-3`)]

hetLsT <- lapply(hetLs,function(x){strsplit(x,split = ":")})
hetLsT <- lapply(hetLsT,function(x){(lapply(x,function(y){y[2]}))})
hetLsT <- lapply(hetLsT,function(x){lapply(x,strsplit,split = ",")})
hetRatLs <- lapply(hetLsT,function(x){unlist(lapply(x,function(y){y<-as.numeric(unlist(y)[1:2]);min(y)/sum(y)}))})
hetSumLs <- lapply(hetLsT,function(x){unlist(lapply(x,function(y){y<-as.numeric(unlist(y)[1:2]);sum(y)}))})

plot(hetSumLs[[1]],hetRatLs[[1]])

hist(hetRatLs[[1]],1000)
hist(hetRatLs[[2]],1000)
hist(hetRatLs[[3]],1000)
hist(hetRatLs[[4]],1000)

hist(table(sample(c(1:100000,1:100000,1:100000),size = 100000*25.5,replace = T)),1000)


