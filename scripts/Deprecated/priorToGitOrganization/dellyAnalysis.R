#Setup environment
setwd("/home/SharedUser/FastNeutron")
vcfFile<- "Plate1-9.sc2.delly.merged.filtered.anno.vcf"
outBaseName<-"delly"

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

#Format full header
fullHeader<-strsplit(tail(vcfHeader,n = 1),split = "\\t")[[1]]

###View(cbind(fullHeader,longestInfoRow.split))

#Full data analysis
vcfFileLength<-as.numeric(gsub(" .*","",system(paste0("wc -l ",vcfFile),intern = T)))
system(paste0("tail -n ",vcfFileLength-tableBeg+1," ",vcfFile," > ",outBaseName,"fullTable.tsv"),wait = T)

fullVcfTable <- read.table(paste0(outBaseName,"fullTable.tsv"),stringsAsFactors = F,comment.char = "",header= T)

#isolate unique vars
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
splitInfo    <-strsplit(fullVcfTable$INFO,fixed = T,split = ";")
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

#Reassemble
fullTableSummary <- cbind(vcfCut2,splitInfo.df)
fullTableSummary$infoField1<-gsub(";.*","",fullTableSummary$INFO)
fullTableSummary$ncharCons<-nchar(fullTableSummary$CONSENSUS)
fullTableSummary$EndMinusPoSMinus1<-as.numeric(fullTableSummary$END)-as.numeric(fullTableSummary$POS)-1

#Hist
delMutLengthsInBp<-fullTableSummary$EndMinusPoSMinus1[fullTableSummary$SVTYPE=="DEL"]
hist(log10(delMutLengthsInBp),main=paste0(length(delMutLengthsInBp)," del typed rows"),breaks = 1000)

# table output
write.csv(fullTableSummary,file = paste0(outBaseName,"FullSummary.csv"))






