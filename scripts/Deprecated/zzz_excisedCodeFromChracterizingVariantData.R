
#Manually set known relationships based on best-judgement
table(samDf$father,useNA = "ifany")


uniFatherless <- samDf$Group.1[is.na(samDf$father)]
df1Homo <- df1[df1$var_count==2,]
for(i in 1:length(uniFatherless)){
  if(i==1){
    basDf <- NULL
  }
  currVars  <- df1[df1$sample==uniFatherless[i],]
  currAgg     <- aggregate(df1$id    ,by=list(df1$sample    ),function(x,y){sum(x%in%y)},y=currVars$id)
  currHomoAgg <- aggregate(df1Homo$id,by=list(df1Homo$sample),function(x,y){sum(x%in%y)},y=currVars$id)
  currAgg$homo_x <- 0
  currAgg$homo_x [match(currHomoAgg$Group.1,currAgg$Group.1)]<- currHomoAgg$x
  currAgg$fatherlessName   <- currVars$sample[1]
  basDf<- rbind(basDf,currAgg)
  print(i)
}
basDf <- basDf[basDf$Group.1!=basDf$fatherlessName,]
basDf <- basDf[basDf$x>0,]

basDf$flMother <- df1$mother[match(basDf$fatherlessName,df1$sample)]
basDf$pfMother <- df1$mother[match(basDf$Group.1       ,df1$sample)]
basDf$pfFather <- samDf$father[match(basDf$Group.1       ,samDf$Group.1)]
basDf$n12      <- samDf$lfm_n[match(basDf$fatherlessName,samDf$Group.1)]
basDf$n11      <- samDf$lfm_n[match(basDf$Group.1       ,samDf$Group.1)]

basDf$inbredToSib <- (basDf$n11*0.5-basDf$x)^2/(basDf$n11*0.5)
hist(basDf$inbredToSib,1000)


basDf$n1 <- basDf$n11/0.75
basDf$n2 <- (basDf$n12 - 0.5*(basDf$n1))/0.5

basDf$n2 - basDf$n1

# P1xP2  vs. P1xP1 assuming P1 fully heterozygous and P1 and P2 have similar frequencies of mutation
#P1xP2 heterozygous at half of P1 
#P1xP2 heterozygous at half of P1 

#P1xP1 homozygous mutant at 33% of P1
#P1xP1 heterozygous at 

View(basDf[order(-basDf$x),])

basDf$x

#Look for 
df1VarTable <-aggregate(df1$var_count,by=list(df1$sample),table)
df1VarFreqs <- do.call("rbind",lapply(df1VarTable$x,function(x){c(het=x["1"],homo=x["2"])}))
df1VarFreqs[is.na(df1VarFreqs)]<-0
plot(df1VarFreqs[,1],df1VarFreqs[,2],col=c("black","red")[(df1VarTable$Group.1%in%(df1$sample[df1$outcross]))+1])

heatmap(samMat[samplesInvestigated,samplesInvestigated],scale = "none")

#Explore the results
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
