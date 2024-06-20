#dir.create("~/Experiments/PhalFNP/data_ignored/secondary/",recursive = T)
setwd("~/Experiments/PhalFNP/")

library(data.table)

df1 <- data.table::fread("data_ignored/secondary/top1Mlines_depth.tsv",header = F)

hist(apply(df1,1,median),1000)

#Initial position removal
mode<-function(x){
    as.numeric(names(which.max(table(x))))
}
df1_mode <- apply(df1[,3:ncol(df1)],2,mode)

df1_modeCentered <- cbind(df1[,1:2],t(t(df1[,3:ncol(df1)]) - df1_mode))

x <- unlist(df1[,3])
halfPeakTrans <- function(x,center=T){
    x_dens <- density(x,bw = "SJ",kernel = "triangular")
    x_dens$y[x_dens$x==min(x_dens$x)]<-0
    x_peak_ind    <- which.max(x_dens$y) 
    x_peak_height <- x_dens$y[x_peak_ind]
    x_peak_val    <- x_dens$x[x_peak_ind]
    
    x_left_logic  <- x_dens$x<x_peak_val
    x_right_logic <- x_dens$x>x_peak_val
    
    x_leftHalfPeak  <- x_dens$x[x_left_logic ][which.min(abs(x_dens$y[x_left_logic ] - 0.5*x_peak_height))]
    x_rightHalfPeak <- x_dens$x[x_right_logic][which.min(abs(x_dens$y[x_right_logic] - 0.5*x_peak_height))]
    
    x_rightMinusLeft <- x_rightHalfPeak - x_leftHalfPeak
    x_MeanLeftRight <- mean(c(x_leftHalfPeak, x_rightHalfPeak))
    
    #plot(x_dens,xlim=c(x_leftHalfPeak-20,x_rightHalfPeak+20))
    #abline(h = 0.5*x_peak_height,col="red")
    #abline(v = c(x_leftHalfPeak,x_MeanLeftRight,x_rightHalfPeak),col="blue")
    
    if(center){
        return((x - x_MeanLeftRight)/x_rightMinusLeft)
    }else{
        return(x/x_MeanLeftRight)
    }
}

#Pass 1
filterFunction <- function(dfF,keptQuant = 0.5){
    dfC <- apply(dfF[,3:ncol(dfF)],2,halfPeakTrans)
    rAbsMeans <- abs(apply(dfC,1,median))
    #hist(abs(rowMeans1),2000,xlim=c(-2,2))
    cutOff <- quantile(abs(rAbsMeans),probs = keptQuant)
    #abline(v = quantile(abs(rowMeans1),probs = 0.01),col = "red")
    return(dfF[rAbsMeans<cutOff,])
}
df2 <- filterFunction(df1) #Half
df3 <- filterFunction(df2) #1/4
df4 <- filterFunction(df3) #1/8
df5 <- filterFunction(df4) #1/16

plot(x = 1:max(df5$V2),ylim=c(0,6),type="n")
df5_trans <- apply(df5[,3:ncol(df5)],2,halfPeakTrans,center = F)
colVec <- rainbow(ncol(df5_trans))
colVec <- sample(colVec)
i<-1
for(i in 1:ncol(df5_trans)){
    points(df5$V2,df5_trans[,i],col=colVec[i])
    print(i)
}
