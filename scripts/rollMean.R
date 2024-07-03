wd <- "~/Experiments/PhalFNP/"
prelimFile <- "data_ignored/secondary/prelimAggregateDepths_20240625_1442.csv"

#Load data
setwd(wd)
library(data.table)
df1 <- fread("data_ignored/secondary/aggregateDepths4.csv")

library(zoo)

df1 <- df1[order(df1$ind,df1$chr,df1$minInt),]
df1$ind_chr <- paste0(df1$ind,"_",df1$chr)
setkey(df1,ind_chr)

dtAgg <- df1[, frollmean(avg_FullNorm_Trimmed, n = 3,align="center" adaptive=T, by = ind_chr)]


df1$prevNorm <- c(NA,df1$avg_FullNorm_Trimmed[1:(nrow(df1)-1)])
df1$nextNorm <- c(df1$avg_FullNorm_Trimmed[2:(nrow(df1))],NA)

chrLogicPrev <- c(df1$chr[2:nrow(df1)]!=df1$chr[1:(nrow(df1)-1)],F)
chrLogicNext <- c(F,df1$chr[2:nrow(df1)]!=df1$chr[1:(nrow(df1)-1)])
df1$prevNorm[chrLogicPrev]<-NA
df1$nextNorm[chrLogicNext]<-NA

df1$win3 <- (df1$nextNorm +df1$avg_FullNorm_Trimmed + df1$prevNorm)/3

#View(df1[1465:1485,])

df1Plotted <- df1[grepl("Chr06",df1$chr),]
library
p1 <- ggplot(df1Plotted,
             aes((minInt+maxInt)/2,win3,color=ind,group=ind))+
  geom_line()+
  theme_bw()+
  facet_wrap(~chr,ncol = 2)+
  guides(color="none",group="none")+
  scale_y_continuous(limits = c(-15,15))
p1
