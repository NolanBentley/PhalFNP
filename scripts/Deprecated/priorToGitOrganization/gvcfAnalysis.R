#Setup environment
setwd("/home/SharedUser/FastNeutron")
vcfFile<- "Plate1-9.sc2.genotype_gvcfs.f1.bf=g10-G3-Q40-QD5-MQ40-GQ20.anno.vcf"
outBaseName<-"gvcfs"

#Use a bash command to operate on the file stream extremely efficiently and "cut" the data
## Note that the string in the following is actually a line of code in a different computing language called bash
## This won't work on a non bash using system (like most windows machines)
system(paste0("cut -f 1,2,3,4,5,6,7,8,9,10 ",vcfFile," > ",outBaseName,"Cut1-10.vcf"),wait = T)

#Read the subset data into R
## Don't try to print! It ias too big
vcfCut <- readLines(paste0(outBaseName,"Cut1-10.vcf"))

#Remove header
tableBeg  <- grep("#CHROM",vcfCut)[1]
vcfHeader <- readLines(vcfFile,n = tableBeg)
vcfCut    <- vcfCut[tableBeg:length(vcfCut)]

#Print table header 
vcfCut[1]
vcfCut[1]<- gsub("^#","",vcfCut[1])#Remove # to prevent comment char interpretation

#Save file
write(vcfCut,file = paste0(outBaseName,"Cut1-12.tsv"))
#Read table
vcfCut2<-read.table(paste0(outBaseName,"Cut1-12.tsv"),header = T,stringsAsFactors = F)
vcfCut2

head(vcfCut2)

hist(nchar(vcfCut2$INFO),breaks = 1000)

precise<- vcfCut2$INFO[nchar(vcfCut2$INFO)>140&nchar(vcfCut2$INFO)<200]

longestInfoRow <- which.max(nchar(vcfCut2$INFO)) + tableBeg
longestInfoRow.raw<- system(paste0("head -",longestInfoRow," ",vcfFile," | tail -n 1"),intern = T)
longestInfoRow.split<-as.character(strsplit(longestInfoRow.raw,split = "\\t")[[1]])

#Check my efforts
longestInfoRow.split[8]==vcfCut2$INFO[which.max(nchar(vcfCut2$INFO))]

#Format full header
fullHeader<-strsplit(tail(vcfHeader,n = 1),split = "\\t")[[1]]

###View(cbind(fullHeader,longestInfoRow.split))

#Full data analysis
vcfFileLength<-as.numeric(gsub(" .*","",system(paste0("wc -l ",vcfFile),intern = T)))
system(paste0("tail -n ",vcfFileLength-tableBeg+1," ",vcfFile," > ",outBaseName,"fullTable.tsv"),wait = T)

fullVcfTable <- read.table(paste0(outBaseName,"fullTable.tsv"),stringsAsFactors = F,comment.char = "",header= T)

#Check for no repeats
fullTable.strings<-paste0(fullVcfTable$X.CHROM,fullVcfTable$POS,fullVcfTable$ID)
if(!max(table(fullTable.strings))==1){stop("Oh no!")
}else{rownames(fullVcfTable)<-fullTable.strings}

#isolate uinique vars
varTable<-fullVcfTable[,10:ncol(fullVcfTable)]
varTable <- apply(varTable,2,gsub,pattern="\\:.*",replacement="")

isMut<-apply(varTable,2,`%in%`,c("0/1","1/1"))
fullVcfTable$numGenos<-rowSums(isMut)
fullVcfTable$lines<-NA
for(i in 1:nrow(isMut)){
  fullVcfTable$lines[i]<-paste0(colnames(isMut)[isMut[i,]],":",varTable[i,][isMut[i,]],collapse=",")
  vcfCut2$Lines[i]<-fullVcfTable$lines[i]
  if(i%%1000==1){print(i)}
}

#Visualize mutation frequency
hist(colSums(isMut),100)

#Split the mutation info
splitInfo<-strsplit(fullVcfTable$INFO,fixed = T,split = ";")
splitInfo.cat<-lapply(splitInfo,gsub,pattern = "=.*",replacement="")
splitInfo.val<-lapply(splitInfo,gsub,pattern = ".*=",replacement="")

#Create a table of all category values
catFreq<-table(unlist(splitInfo.cat))
splitInfo.df<-data.frame(matrix(NA,nrow = nrow(fullVcfTable),ncol = length(catFreq)))
rownames(splitInfo.df)<-rownames(fullVcfTable)
colnames(splitInfo.df)<-names(catFreq)[order(-catFreq)]
for(i in 1:nrow(splitInfo.df)){
  matchingColumns<-match(splitInfo.cat[[i]],colnames(splitInfo.df))
  splitInfo.df[i,matchingColumns]<-splitInfo.val[[i]]
  if(i%%1000==1){print(i)}
}

splitInfo.df$Impact<-gsub("\\(.*","",splitInfo.df$EFF)

cbind(names(sort(colSums(isMut[splitInfo.df$Impact!="INTERGENIC",])))[1:20],
            sort(colSums(isMut[splitInfo.df$Impact!="INTERGENIC",])) [1:20])

#Reassemble
fullTableSummary <- cbind(vcfCut2,splitInfo.df)
fullTableSummary$infoField1<-gsub(";.*","",fullTableSummary$INFO)

# table output
write.csv(fullTableSummary,file = paste0(outBaseName,"FullSummary.csv"))

###View(fullTableSummary[fullTableSummary$SVTYPE%in%c("DEL","INS")&fullTableSummary$infoField1=="PRECISE",])


