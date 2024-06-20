#Setup environment
setwd("/home/SharedUser/FastNeutron")

#Download the gff3 file
## THis wont work again.
#system("curl --cookie jgi_session=/api/sessions/30b40ffc0893003711e9c1a2f3c634ae --output download.20220616.150633.zip -d \"{\"ids\":{\"Phytozome-495\":[\"59cd36907ded5e2f18695aeb\"]}}\" -H \"Content-Type: application/json\" https://files.jgi.doe.gov/filedownload/",intern = T)
#system("unzip download.20220616.150633.zip",intern = T)
#system("gunzip Phytozome/PhytozomeV12_unrestricted/Phallii/annotation/Phallii_495_v3.1.gene_exons.gff3.gz",intern = T)
#make this cp instead #system("mv Phytozome/PhytozomeV12_unrestricted/Phallii/annotation/Phallii_495_v3.1.gene_exons.gff3 ./",intern = T)

#Read gff3 file
annoFile<-"Phallii_495_v3.1.gene_exons.gff3"
phAnno <- read.table(file = annoFile,header = F,fill = T,stringsAsFactors = F)
phAnno$rowNumber<-1:nrow(phAnno)

#What is up with ;?
#phAnnoWeird<-which(phAnno$V3==";")
#View(phAnno[sort(unique(c(phAnnoWeird-1,phAnnoWeird,phAnnoWeird+1))),])
##It is fine...

phAnno.exons<-phAnno#[phAnno$V3=="exon",]
phAnno.exons$geneId<-gsub("...v3\\.1.*","",gsub("ID\\=","",gsub("\\;.*","",phAnno.exons$V9)))
phAnno.exons$mutsWithin<-NA


#Read a vcf file summary
vcf.file <- "dellyFullSummary.csv"
vcfSum <- read.csv(vcf.file,stringsAsFactors = F)

uniChr <- unique(vcfSum$CHROM)

vcfSum$withinGenes <- NA
vcfSum$rowNumbers  <- NA

currChrInd <- 1
currRow <- 1
for(currChrInd in 1:length(uniChr)){
  currChr<-uniChr[currChrInd]
  currVcf <-vcfSum[vcfSum$CHROM==currChr,]
  currAnno<-phAnno.exons[phAnno.exons$V1==currChr,]
  
  #Put genes in muts
  for(currRow in 1:nrow(currVcf)){
    currExtents<-c(currVcf$POS[currRow],currVcf$END[currRow])
    currGenes<-which(!(max(currExtents,na.rm=T)<currAnno$V4|min(currExtents,na.rm=T)>currAnno$V5))
    currVcf$withinGenes[currRow]<-paste0(unique(currAnno$V9[currGenes]),collapse=";")
    currVcf$rowNumbers [currRow]<-paste0(currAnno$rowNumber[currGenes],collapse = ";")
  }
  
  #Subset searched for genes
  genesSearches <- sort(as.numeric(unique(unlist(strsplit(unique(currVcf$rowNumbers),split = ";")))))
  rowList       <- lapply(strsplit(unique(currVcf$rowNumbers),split = ";"),as.numeric)
  
  #Put mut lines in genes
  if(length(genesSearches)>0){
    for(currRow in 1:length(genesSearches)){
      currMuts<-which(unlist(lapply(rowList,function(x){any(genesSearches[currRow]%in%x)})))
      currExonRow<-match(genesSearches[currRow],currAnno$rowNumber)
      currAnno$mutsWithin[currExonRow]<-paste0(unique(currVcf$Lines[currMuts]),collapse=";")
      if(currRow%in%round(seq(1,length(genesSearches),length.out = 100))){cat(".")}
    }
  }else{cat("|")}
  
  #Save
  if(length(genesSearches)>0){
    print(currChrInd)
    phAnno.exons[phAnno.exons$V1==currChr,]<-currAnno
  }
  vcfSum[vcfSum$CHROM==currChr,]<-currVcf
}

write.csv(vcfSum,paste0(gsub("\\.csv$","",vcf.file),"_WithAnnotation.csv"))
write.csv(phAnno.exons,paste0(gsub("\\.gff3$","",annoFile),"_WithMutations.csv"))

#Analyze
phAnno.exons<-read.csv(paste0(gsub("\\.gff3$","",annoFile),"_WithMutations.csv"),stringsAsFactors = F)

library(openxlsx)#install.packages("openxlsx")
openxlsx::write.xlsx(phAnno.exons, paste0(gsub("\\.gff3$","",annoFile),"_WithMutations.xlsx"), sheetName = "Sheet1")

#Cluster the samples
anno.mut<-phAnno.exons[!is.na(phAnno.exons$mutsWithin),]
anno.mut$trimmedMuts<-gsub("\\:.\\/.","",anno.mut$mutsWithin)
mutList<-strsplit(anno.mut$trimmedMuts,split = ",|;")

length(mutList)

uniMut<-sort(unique(unlist(mutList)))

anno.mat<-matrix(0,length(anno.mut$geneId),length(uniMut),
                 dimnames = list(anno.mut$geneId,uniMut))
i<-1
for(i in 1:nrow(anno.mat)){
  if(length(mutList[[i]])>0){
    anno.mat[i,mutList[[i]]]<-1
    stop()
  }
  if(i%%1000==1){print(i)}
}

hist(rowMeans(anno.mat),1111)
geneMutFreq<-rowSums(anno.mat)
geno.dist<-dist(t(anno.mat[geneMutFreq>3,]))
save(geno.dist,file = "geno.dist.rdata")
geno.cmd <-cmdscale(geno.dist)

#install.packages("heatmaply")
require(heatmaply)


plot(geno.cmd[,1],geno.cmd[,2])



