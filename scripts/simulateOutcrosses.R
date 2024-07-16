## Lead functions
source("~/Experiments/PhalFNP/scripts/similarityCalc.R")

## Enable for testing
combineWithDf1 <- T

#Simulate progeny
nPar     <- 7
minVars  <- 10
multBtwn <- 20
reps     <- 5
for(j in 1:nPar){
  if(j==1){
    dfSim <- NULL
    currRows<-0
  }else{
    currRows <- nrow(dfSim)
  }
  dfSim <- rbind(dfSim,data.frame(
    sample=paste0("simPar",j),
    id=paste0("simVar",(currRows+1):(currRows+minVars+multBtwn*(j-1))),
    var_count=1
  ))
}
table(dfSim$sample)
combos <- expand.grid(unique(dfSim$sample),unique(dfSim$sample))
combos$prog <- paste0("simProg",gsub("simPar","",combos$Var1),"x",gsub("simPar","",combos$Var2))
for(k in 1:nrow(combos)){
  for(i in 1:reps){
    dfSim <- rbind(dfSim,makeProgeny(dfSim,combos$Var1[k],combos$Var2[k],paste0(combos$prog[k],"-",i)))
  }
}
dfSim <- dfSim[dfSim$var_count>0,]



length(unique(dfSim$sample))
dfSimVarTable <-aggregate(dfSim$var_count,by=list(dfSim$sample),table)
lapply(dfSimVarTable$x,function(x){c(het=x["1"],homo=x["2"])/sum(x)})

#Merge and find duplicates for parsing
if(combineWithDf1){
  df2 <- merge(df1,dfSim,all = T)
}else{
  df2 <- dfSim
}

#Calculate values for loop
relaDf <- similarityFun(df2)

#Make long form similarity
relaDf$logicSum <-rowSums(relaDf[,grep("^(het|homo)_(het|homo)$",colnames(relaDf))])
relaDf$dupeSum  <-rowSums(relaDf[,grep("len_(i|j)",colnames(relaDf))])
relaDf$totSum   <-rowSums(relaDf[,grep("tot_(i|j)",colnames(relaDf))])
relaDf$totProp  <- relaDf$logicSum/relaDf$totSum

#Make smaller
hist(relaDf$totProp,ylim=c(0,1000),10000);abline(v = 0.025,col="red")
relaDf <- relaDf[relaDf$totProp>0.025,]

#Add in columns at beginning
relaDf <- cbind(comp=trimws(mapply(function(x,y){return(paste0(sort(c(x,y)),collapse = "|"))},format(relaDf$iSam,justify="left"),format(relaDf$jSam,justify="left"))),relaDf[,!grepl("^comp",colnames(relaDf))])
relaDf <- cbind(simRela="Unknown",relaDf[,!grepl("simRela",colnames(relaDf))])

#Filter and organize
relaDf <- relaDf[!duplicated(relaDf$comp),]
relaDf <- relaDf[relaDf$iSam!=relaDf$jSam,]
relaDf <- relaDf[order(-relaDf$totProp),]

#View and add in known relationships
relaDf$comp_ped <- gsub(" ","",gsub("-[[:digit:]]+","",relaDf$comp))
relaDf$bothSim <- grepl("^sim",relaDf$iSam)&grepl("^sim",relaDf$jSam)
cat(paste0('relaDf$simRela[relaDf$comp_ped=="',unique(relaDf$comp_ped[relaDf$simRela=="Unknown"&relaDf$bothSim]),'"]<-""',collapse = "\n"))
relaDf$areIdentical <- relaDf$iSam==relaDf$jSam
relaDf$iGeno <- gsub("simProg|-.*","",gsub("simPar","0x",relaDf$iSam))
relaDf$jGeno <- gsub("simProg|-.*","",gsub("simPar","0x",relaDf$jSam))
relaDf$iA    <- substr(relaDf$iGeno,1,1)
relaDf$iB    <- substr(relaDf$iGeno,3,3)
relaDf$jA    <- substr(relaDf$jGeno,1,1)
relaDf$jB    <- substr(relaDf$jGeno,3,3)
relaDf$iPar  <- grepl("Par",relaDf$iSam)
relaDf$jPar  <- grepl("Par",relaDf$jSam)
relaDf$iSelf <- relaDf$iA==relaDf$iB
relaDf$jSelf <- relaDf$jA==relaDf$jB
relaDf$ijSum <- mapply(function(x,y,a,b){sum(c(x,y)%in%c(a,b))},relaDf$iA,relaDf$iB,relaDf$jA,relaDf$jB)
relaDf$jiSum <- mapply(function(x,y,a,b){sum(c(x,y)%in%c(a,b))},relaDf$jA,relaDf$jB,relaDf$iA,relaDf$iB)
relaDf$ijDescriptor <- paste(relaDf$areIdentical,relaDf$iPar,relaDf$jPar,relaDf$iSelf,relaDf$jSelf,relaDf$ijSum,relaDf$jiSum,sep = ".")
relaDf$ijDescriptor[!grepl("sim",relaDf$iSam)]<-"Unknown"

descTable <- table(relaDf$ijDescriptor[relaDf$simRela=="Unknown"])
cat(paste0("relaDf$comp[relaDf$ijDescriptor=='",names(descTable),"'] # ",descTable,collapse = "\n"))
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.FALSE.1.1'] # 92849
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.FALSE.2.2'] # 9135
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.TRUE.1.2'] # 9450
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.FALSE.2.1'] # 9447
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.TRUE.2.2'] # 735
relaDf$comp[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.FALSE.1.1'] # 1260
relaDf$comp[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.TRUE.1.2'] # 105

relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.FALSE.0.0'] # 240
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.TRUE.0.0'] # 120
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.FALSE.0.0'] # 120
relaDf$comp[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.TRUE.0.0'] # 40
relaDf$comp[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.FALSE.0.0'] # 120
relaDf$comp[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.TRUE.0.0'] # 40
relaDf$comp[relaDf$ijDescriptor=='FALSE.TRUE.TRUE.FALSE.FALSE.1.1'] # 10

descTable <- table(relaDf$ijDescriptor[relaDf$simRela=="Unknown"])
cat(paste0("relaDf$simRela[relaDf$ijDescriptor=='",names(descTable),"'] <- '' # ",descTable,collapse = "\n"))
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.FALSE.1.1'] <- 'Half-siblings' # 92849
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.FALSE.2.2'] <- 'Full-siblings' # 9135
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.TRUE.1.2' ] <- 'Sibling to Self' # 9450
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.FALSE.2.1' ] <- 'Self to sibling' # 9447
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.TRUE.2.2'  ] <- 'Self to self' # 735
relaDf$simRela[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.FALSE.1.1' ] <- 'Parent / child' # 1260
relaDf$simRela[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.TRUE.1.2'  ] <- 'Parent / self' # 105
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.FALSE.0.0'] <- 'UR children' # 240
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.FALSE.TRUE.0.0'] <- 'UR child/self' # 120
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.FALSE.0.0'] <- 'UR self/child' # 120
relaDf$simRela[relaDf$ijDescriptor=='FALSE.FALSE.FALSE.TRUE.TRUE.0.0'] <- 'UR selfs' # 40
relaDf$simRela[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.FALSE.0.0'] <- 'Parent / UR cross' # 120
relaDf$simRela[relaDf$ijDescriptor=='FALSE.TRUE.FALSE.FALSE.TRUE.0.0'] <- 'Parent / UR self' # 40
relaDf$simRela[relaDf$ijDescriptor=='FALSE.TRUE.TRUE.FALSE.FALSE.1.1'] <- 'Parent / Parent' # 10
table(relaDf$simRela)

#Remove parents since that should be possible
relaDf <- relaDf[!grepl("Parent",relaDf$simRela),]
#relaDf <- relaDf[order(relaDf$simRela=="Unknown",relaDf$totProp),]

#Add in proportional analyses
relaDf$minTot  <- apply(cbind(relaDf$tot_i,relaDf$tot_j),1,min)
relaDf$meanTot <- (relaDf$tot_i+relaDf$tot_j)/2
propDf <- relaDf[,grep("^(het|homo)_(het|homo)$",colnames(relaDf))]
colnames(propDf) <- paste0(colnames(propDf),"_trans")
relaDf <- cbind(relaDf[,!grepl("_trans$",colnames(relaDf))],propDf)

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
relaDf$crToIb_chiSum <- chiFromTot(relaDf,relaDf$tot_j,c(3,6,0,0))
relaDf$ibToCr_chiSum <- chiFromTot(relaDf,relaDf$tot_i,c(3,0,6,0))
relaDf$ibToIb_chiSum <- chiFromTot(relaDf,(relaDf$tot_j+relaDf$tot_i)/2,c(3,6,6,12))
relaDf$fs_chiSum <- chiFromTot(relaDf,relaDf$tot_i,c(2,0,0,0))+
  chiFromTot(relaDf,relaDf$tot_j,c(2,0,0,0))

if(combineWithDf1){
  write.csv(relaDf,file = "data/obsAndSimRelas.csv")
}else{
  write.csv(relaDf,file = "data/simulatedRelas.csv")
}


library(ggplot2)
hist(relaDf$totProp,1000)
ggplot(relaDf[relaDf$totProp>0.025,],aes(inbredToSib,totProp,color=simRela,shape=simRela))+
  geom_point()


relaDf$het_het-relaDf$totSum*4

#Save analysis
save.image(file = "./data_ignored/secondary/SimulationPriorToVis.rimage")

#PCA
library(ggplot2) 
library(ggfortify)
pcaData_withGrp <- relaDf[,grep("_trans|^simRela$",colnames(relaDf))]
data.frame(
  fullSib <- relaDf$het_het
)
pcaData <- log10(pcaData_withGrp[,colnames(pcaData_withGrp)!="simRela"]+1)/log10(relaDf$totSum+1)
pca<-prcomp(pcaData,scale. = T)
autoplot(pca,data=pcaData_withGrp,
         color="simRela",
         shape="simRela",
         loadings=T,loadings.label=T)+
  theme_bw()




i<-1
numCols <- which(lapply(relaDf,class)=="numeric")
for(i in 1:length(numCols)){
  currAgg <- aggregate(relaDf[,numCols[i]],by=list(relaDf$simRela),mean)
  if(i == 1){
    sumAgg <- data.frame(group=currAgg$Group.1,newCol=currAgg$x)
  }else{
    sumAgg$newCol<-currAgg$x
  }
  colnames(sumAgg)[i+1]<-colnames(relaDf)[i]
}
sumAgg

save.image("reladDf.rimage")
bplottedRela <-relaDf
rownames(plottedRela) <- relaDf$plottedNames
plottedRela <- plottedRela[order(plottedRela$plottedNames),sapply(plottedRela,class)=="numeric"]

library(heatmaply)
library(htmlwidgets)
hm1 <- heatmaply(as.matrix(plottedRela),Rowv = F,Colv = F)
saveWidget(hm1,file = "data/simHm.html")

#Save similarity data
write.csv(samMat,file = "~/Experiments/PhalFNP/data/similarityMatrix.csv")
write.csv(samMat,file = "~/Experiments/PhalFNP/data/similarityLong.csv")