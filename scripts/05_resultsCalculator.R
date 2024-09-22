#Setup
wd <- "~/Experiments/PhalFNP/"; setwd(wd)
coverageFile <- "data/sampleDf.csv"
variantFile  <- "https://utexas.box.com/shared/static/yht8ojfafsab3btuq6xe2chw8by9onfo.csv"

#Load data
samDf <- read.csv(coverageFile)
df1   <- read.csv(variantFile )

#Mutation burden
mean(aggregate(df1$id,by = list(df1$sample),length)$x)
