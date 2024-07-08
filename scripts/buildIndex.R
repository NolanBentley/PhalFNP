#### Load aggDf ####
wd <- "~/Experiments/PhalFNP/"
setwd(wd)
if(!exists("aggDf")){aggDf<-read.csv("data_ignored/secondary/plottedAggDf_n.csv")}

#### Build index ####
index <- c(
  '<!DOCTYPE html>','<html>','<head>','<title>',
  '<i>Panicum hallii</i> fast-neutron population large structural variant explorer',
  '</title>','</head>','<body>','<h1>',
  'Welcome to the large variant explorer!',
  '</h1>','<p>',
  'This website houses interactive plots for exploring coverage-based variants from the P. hal FNP.',
  
  '<br><br>As of 7/4/2024, this data visualizes two primary statistics; <br>',
  '<b>Hexagons:</b> the frequency of bins of position and normalized-coverage across the population being visualized.<br>',
  '<b>Points:</b> the normalized-coverage for specific samples subset to only the intervals where the absolute value of the y-axis exceeds certain cutoffs.',
  
  '<br><br>The y-axis is calculated from the median coverage transformed into estimated ploidy across 5 intervals. ',
  'Each interval summarizes the mean read depth across a 2.5-5kb range of positions excluding positions with 0 depth. ',
  '<br>The coverage has been transformed to estimate copy-number by multiplying by 2 and dividing by the peak value.',
  '</p>','</body>','</html>'
)
index_endBody <- which(trimws(index)=="</body>")
aggDf$singleLineBinPages <- paste0("./depthImages/singleLineIndexes/",aggDf$idNum_bin,".html")
singleLineBinPages <-sort(unique(aggDf$singleLineBinPages))
index_added <- c(
  "<br><br><b><u>multi-sample multi-chromosome plots</b></u>",
  "<br>Clicking on points in these plots will take you to that chromosome's multi-sample plot.",
  "<br>Note: These are large (~40Mb) files and will take a bit to load and render.",
  "<ul>",
  sort(unique(aggDf$multiLineMultiChr)),
  "</ul><br><br><b><u>multi-sample single-chromosome plots</b></u>",
  "<br>Clicking on points in these plots will take you to that sample's multi-chromosome plot. Clicking on those points will take you to the single-chromosome plot for that line.",
  "<ul>",
  sort(unique(aggDf$multiLineSingleChr)),
  "</ul><br><br><b><u>Single line plots</b></u>",
  "This section takes you to single-sample plots based on the sample ID.",
  "<ul>",
  singleLineBinPages,"</ul>"
)
htmlPos <-grep("html",index_added)
index_mod <- index_added
index_mod[htmlPos] <- paste0('<li><a href="',
                             gsub("^\\.","/PhalFNP",index_added[htmlPos]),
                             '">',basename(index_added[htmlPos]),
                             "</a></li>"
)
index_built <- c(
  index[1:(index_endBody-1)],
  index_mod,
  index[index_endBody:length(index)]
)
write(index_built,"index.html")

#Make the position html files
posHtmlFiles <- unique(index_added[htmlPos])
posHtmlFiles <- posHtmlFiles[grep("Lines",posHtmlFiles)]
dir.create(dirname(posHtmlFiles[1]),showWarnings = F)
i<-1
for(i in 1:length(posHtmlFiles)){
  index <- c(
    '<!DOCTYPE html>','<html>','<head>','<title>',
    basename(posHtmlFiles),
    '</title>','</head>','<body>','<h1>',
    'Welcome to the large variant explorer for ',basename(posHtmlFiles),
    '</h1>','<p>',
    '</p>','</body>','</html>'
  )
  index_endBody <- which(trimws(index)=="</body>")
  currRows <- aggDf$singleLineBinPages == posHtmlFiles[i]
  index_added <- c(
    "<br><br><b><u>single-sample multi-chromosome plots</b></u>",
    "<br>Clicking on points in these plots will take you to that chromosome's single-sample plot.",
    "<ul>",
    sort(unique(aggDf$singleLineMultiChr[currRows])),
    "</ul><br><br><b><u>single-sample single-chromosome plots</b></u>",
    "<br>Clicking on points in these plots will take you to that sample's multi-chromosome plot. Clicking on those points will take you to the single-chromosome plot for that line.",
    "<ul>",
    sort(unique(aggDf$singleLineSingleChr[currRows]))
  )
  htmlPos <-grep("html",index_added)
  index_mod <- index_added
  index_mod[htmlPos] <- paste0('<li><a href="',
                               gsub("^\\.","/PhalFNP",index_added[htmlPos]),
                               '">',basename(index_added[htmlPos]),
                               "</a></li>"
  )
  index_built <- c(
    index[1:(index_endBody-1)],
    index_mod,
    index[index_endBody:length(index)]
  )
  write(index_built,posHtmlFiles[i])
}

#Fix page titles
htmlFiles <- list.files("depthImages/",pattern = "html$",full.names = T,recursive = T)
x <- htmlFiles[1]
replacePlotlyTitles <- function(x){
  currHtml <- readLines(x)
  titlerow<- grep("\\<title\\>plotly\\.*title\\>$",currHtml)
  if(length(titlerow)>0){
    currHtml[titlerow[1]]<-paste0("<title>",basename(x),"</title>")
    write(currHtml,file = x)
  }
  print(x)
}
sapply(htmlFiles,replacePlotlyTitles)

#Add gene link
geneDf <- read.csv("data/geneDf.csv")










