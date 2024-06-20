library(data.table)
library(zoo)
df1<-as.data.frame(fread("~/Experiments/FastNeutron/analysis/andieDepth/andieDepths.out"))
df1_orig<-df1

interestingCols <- c('V105','V139','V176','V312','V380','V495','V712')

#Group the regions
df1_region <- cumsum(c(TRUE,abs(diff(df1$V2))>500))


currData<- unlist(df1[,interestingCols])
hist(currData,1000)
abline(v = mean(currData,trim = 0.4),col="red")
abline(v = mean(currData),col="blue")
plot(density(currData))

find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
     pks <- unlist(pks)
     pks
}


#Remove bad mapping
trimSd <- function(x){
  quants <- quantile(x,probs=c(0.25,0.75))
  x_sort <- sort(x)
  x_trimmed <- x[x>quants[1]|x<quants[2]]
  return(sd(x_trimmed))
}
iqm  <-apply(df1[,-1:-2],MARGIN = 2,FUN = mean,trim=0.25)
q25  <-apply(df1[,-1:-2],MARGIN = 2,FUN = quantile,prob=0.25)
q75  <-apply(df1[,-1:-2],MARGIN = 2,FUN = quantile,prob=0.75)
df1_iqmCentered <- t((t(df1[,-1:-2])-iqm)/(q75-q25))

iqm_med <- median(rowMeans(df1_iqmCentered))
iqm_mad <- apply(df1_iqmCentered,1,mad)

initRowLogic <- 
  rowMeans(df1_iqmCentered)>(iqm_med-3*iqm_mad)&
  rowMeans(df1_iqmCentered)<(iqm_med+3*iqm_mad)

#Subset based on logic
df2 <- df1[initRowLogic,]
df2_region <- df1_region[initRowLogic]

#Normalize the data
iqm  <-apply(df2[,-1:-2],MARGIN = 2,FUN = mean,trim=0.25)
q25  <-apply(df2[,-1:-2],MARGIN = 2,FUN = quantile,prob=0.25)
q75  <-apply(df2[,-1:-2],MARGIN = 2,FUN = quantile,prob=0.75)
df_normed <- t((t(df2[,-1:-2])-iqm)/(q75-q25))

hist(rowMeans(df1[,-1:-2]),1000)

#heatmap(df_normed,Rowv = NA)

#Loop across normalized data
uniRegions <- unique(df2_region)
for(i in 1:length(uniRegions)){
  if(i==1){out<-NULL}
  currRegion  <- df_normed[df2_region==uniRegions[i],]  
  meanWindows <- rollmean(currRegion[,-1:-2],k = 500)
  out <- rbind(out,meanWindows)
}



