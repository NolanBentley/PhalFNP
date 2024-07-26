setwd("~/Experiments/PhalFNP/")
mac2 <- readRDS("data_ignored/secondary/midFreqMACs.rds")


#Get id in confusing subset
confG <- readLines("data_ignored/primary/confusingGenotypes")
macDf <- as.data.frame(t(mac2))
if(all(confG%in%rownames(macDf))){
  macDf$Confusing <- as.factor((rownames(macDf)%in%confG)*1)
}else{
  stop("ID missmatch")
}

uniLines <- 1:nrow(macDf)#which(macDf$Confusing==1)
relaDf <- matrix(NA,nrow = length(uniLines),ncol=ncol(mac2))
i<-1

mac2Homo <- mac2 == 2
mac2Het  <- mac2 == 1
nRow <- nrow(mac2)
nCol <- ncol(mac2)
j_homo <- colSums(mac2Homo,na.rm = T)
j_het  <- colSums(mac2Het ,na.rm = T)
k <- 1
for(k in 1:length(uniLines)){
  if(k==1){dfOut <- NULL}
  currCol  <- uniLines[k]
  currRows <- which(mac2Homo[,currCol]|mac2Het[,currCol])
  currDf   <- data.frame(
    k        = k,
    i_sam    = colnames(mac2Homo)[currCol],
    j_sam    = colnames(mac2Homo),
    homo_homo= colSums(mac2Homo[currRows,currCol]&mac2Homo[currRows,],na.rm = T),
    homo_het = colSums(mac2Homo[currRows,currCol]&mac2Het [currRows,],na.rm = T),
    het_homo = colSums(mac2Het [currRows,currCol]&mac2Homo[currRows,],na.rm = T),
    het_het  = colSums(mac2Het [currRows,currCol]&mac2Het [currRows,],na.rm = T),
    i_homo   = rep(j_homo[currCol],nCol),
    j_homo   = j_homo,
    i_het    = rep(j_het [currCol],nCol),
    j_het    = j_het
  )
  dfOut <- rbind(dfOut,currDf)
  print(round(k/length(uniLines),3))
}
dfOut$aSum <- dfOut$homo_homo*2 + dfOut$het_homo + dfOut$homo_het + dfOut$het_het*2
dfOut$nSum <- dfOut$i_homo+dfOut$j_homo+dfOut$i_het+dfOut$j_het
dfOut$aPerN <- dfOut$aSum/dfOut$nSum

hist(dfOut$aPerN,ylim=c(0,length(uniLines)),1000)

View(dfOut[which(dfOut$i_sam=="FIL30_85_H_M2"),])

#From demo
set.seed(123)
ind <- sample(2, nrow(iris),
              replace = TRUE,
              prob = c(0.6, 0.4))
training <- iris[ind==1,]
testing <- iris[ind==2,]
linear <- lda(Species~., training)

#Try math
actualLogic <- colnames(mac2)%in%confG
for(i in 1:1000){
  if(i==1){
    permOut <- list()
    currLogic <- actualLogic
  }else{
    set.seed(i)
    currLogic <- sample(actualLogic)
  }
  permOut[[i]] <- 
    rowMeans(mac2[, currLogic],na.rm = T)-
    rowMeans(mac2[,!currLogic],na.rm = T)
  print(i)
}

hist(unlist(permOut[-1]),1000,ylim = c(0,100))
hist(unlist(permOut[1]),1000,border="red",add=T)

#PCA
pca2 <- prcomp(t(mac2[(rowSums(is.na(mac2))==0)&(permOut[[1]]>0),]),scale=F)

plottedDf <- data.frame(pca2$x)
plottedDf$ID <- rownames(plottedDf)
if(all(confG%in%plottedDf$ID)){
  plottedDf$Confusing <- plottedDf$ID%in%confG
}else{
  stop("ID missmatch")
}
library(plotly)
p1 <- ggplot(plottedDf[order(plottedDf$PC3),],
             aes(PC1,PC2,fill=PC3,label=ID,shape=Confusing))+
  geom_point(size=3)+
  #geom_point(data = plottedDf[plottedDf$Confusing,],color="white",alpha=0.5)+
  scale_fill_viridis_c()+
  coord_fixed()+
  theme_bw()
ggplotly(p1)

library(MASS)

