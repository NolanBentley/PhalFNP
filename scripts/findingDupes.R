variantFile <- "~/../Box/Teaching/XPlants/FNPFiles/FNPVariants_20240626.csv"
df1 <- read.csv(variantFile)

dfDupe<-df1[df1$id%in%(df1$id[duplicated(df1$id)]),]

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
