variantFile <- "https://utexas.box.com/shared/static/yht8ojfafsab3btuq6xe2chw8by9onfo.csv"
df1 <- read.csv(variantFile)

df1$analysis <- "GATK"
df1$analysis[is.na(df1$QD)]<-"DELLY"
df1$id      <- paste0(df1$id,":",df1$chrom,"-",df1$pos,"_",df1$ref_allele,"->",df1$alt_allele,"_",df1$analysis)
df1$isDuped <- df1$id%in%(df1$id[duplicated(df1$id)])
dfDupe<-df1[df1$isDuped,]

uniSamples <- sort(unique(dfDupe$sample))
samMat <- matrix(NA,nrow=length(uniSamples),ncol=length(uniSamples))
dimnames(samMat)<-list(uniSamples,uniSamples)
for(i in 1:ncol(samMat)){
    for(j in 1:nrow(samMat)){
        iIds <- unique(dfDupe$id[dfDupe$sample==uniSamples[i]])
        jIds <- unique(dfDupe$id[dfDupe$sample==uniSamples[j]])
        uniIds <- unique(c(iIds,jIds))
        samMat[i,j]<-(length(iIds)+length(jIds))/2/length(uniIds)
    }
    print(i)
}


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
