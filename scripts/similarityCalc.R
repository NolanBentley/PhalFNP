hasDupes <- function(x){x%in%(x[duplicated(x)])}

cntOfYInX <- function(x,y){
  x_agg <- aggregate(x,by=list(x),length)
  if(!all(x_agg$Group.1%in%x)){stop("Mismatch")}
  z <- rep(NA,length(y))
  z[match(x_agg$Group.1,y)]<-x_agg$x
  return(z)
}

gameteFun <- function(geno){
  geno[geno==1] <- sample(c(0,2),size = sum(geno==1),replace = T)
  return(geno/2)
}
makeProgeny <- function(dfInput, parA, parB, progId){
  parASample <- dfInput[dfInput$sample==parA,]
  parBSample <- dfInput[dfInput$sample==parB,]
  parASample$geno <- gameteFun(parASample$var_count)
  parBSample$geno <- gameteFun(parBSample$var_count)
  prog <- rbind(parASample,parBSample)
  prog_geno <- aggregate(prog$geno,by=list(prog$id),sum)
  prog_out <- data.frame(sample = progId,id=prog_geno$Group.1,var_count=prog_geno$x)
  return(prog_out)
}
iFun <- function(i,samplesInvestigated,uniSamples,dupeList,df2_sample_table,dupeList_var){
  jFun <- function(j,iSam,iIds,iId_cnt,iTot,uniSamples,dupeList,df2_sample_table,dupeList_var){
    jSam <- uniSamples[j]
    jIds <- dupeList[[j]]
    jId_cnt <- dupeList_var[[j]]
    jTot <- df2_sample_table[jSam]
    homo_homo <- jId_cnt==2 & jIds%in%iIds[iId_cnt==2]
    het_homo  <- jId_cnt==2 & jIds%in%iIds[iId_cnt==1]
    homo_het  <- jId_cnt==1 & jIds%in%iIds[iId_cnt==2]
    het_het   <- jId_cnt==1 & jIds%in%iIds[iId_cnt==1]
    
    similarity <- data.frame(
      i=i,j=j,iSam=iSam,jSam=jSam,
      het_het=sum(het_het),
      het_homo=sum(het_homo),
      homo_het=sum(homo_het),
      homo_homo=sum(homo_homo),
      len_i=length(iIds),
      len_j=length(jIds),
      tot_i=iTot,tot_j=jTot
    )
    return(similarity)
  }
  iSam <- uniSamples[i]
  iIds <- dupeList[[i]]
  iId_cnt <- dupeList_var[[i]]
  iTot <- df2_sample_table[iSam]
  out <- do.call("rbind",sapply(samplesInvestigated,FUN = jFun,iSam,iIds,iId_cnt,iTot,uniSamples,dupeList,df2_sample_table,dupeList_var,simplify = F))
  return(out)
}

similarityFun <- function(dfX,ncores=max(c(1,min(c(parallel::detectCores()/2,10))))){
  library(parallel)
  #Calculate values
  dfX_sample_table<-table(dfX$sample)
  dfX$isDuped <- dfX$id%in%(dfX$id[duplicated(dfX$id)])
  dfDupe<-dfX[dfX$isDuped,]
  uniSamples <- sort(unique(dfDupe$sample))
  samplesInvestigated <- 1:length(uniSamples)
  dupeList     <- sapply(uniSamples,function(x){dfDupe$id       [dfDupe$sample==x]})
  dupeList_var <- sapply(uniSamples,function(x){dfDupe$var_count[dfDupe$sample==x]})
  #Loop across data combinations
  cl <- makeCluster(ncores)
  relaDf <- do.call("rbind",parSapply(cl,samplesInvestigated,FUN = iFun,samplesInvestigated,uniSamples,dupeList,dfX_sample_table,dupeList_var,simplify = F))
  stopCluster(cl)
  numCols <- which(colMeans(!is.na(suppressWarnings(apply(relaDf,2,as.numeric))))==1)
  relaDf[,numCols]<- apply(relaDf[,numCols],2,as.numeric)
  return(relaDf)
}





