#Set variables
wd <- "~/Experiments/PhalFNP/"
aggFile <- "data_ignored/secondary/plottedAggDf_n.csv"
toggleFilterToExisting <- T

#### Load aggDf ####
setwd(wd)
if(!exists("aggDf")){aggDf<-read.csv(aggFile)}
aggDf_og <-aggDf
aggDf<-aggDf[aggDf$chrLogic,]

#Remove files if desired
if(toggleFilterToExisting){
  dirsToSearch <- unique(dirname(unique(aggDf$singleLineSingleChr)))
  currFiles <- unique(unlist(sapply(dirsToSearch,list.files)))
  aggDf <- aggDf[basename(aggDf$singleLineSingleChr)%in%currFiles,]
}

#### Build index ####
index <- c('
<!DOCTYPE html>
<html>
<head>
<title>
Panicum hallii fast-neutron population large structural variant explorer
</title>
</head>
<body>
<h1>Welcome to the large variant explorer!</h1>
<p>
This website houses interactive plots for exploring large (>>0.2Mb) copy-number variants from the <i>Panicum. hallii</i> fast neutron population (FNP).<br>
<br>
<dt>The y-axis values were calculated to approximate ploidy across the x-axis intervals:</dt>
<dd>- Each point is the transformed median copy-number across a 5x window of 5Kb intervals every 40Kb for an individual line.</dt>
<dd>- Each interval summarizes the mean read depth across a 2.5-5kb range of positions excluding positions with 0 depth.</dd>
<dd>- The "copy-number" was estimated by multiplying the mean coverage across the interval by 2 and dividing by the peak mean coverage.</dd><br>
<dt>This data is visualized with two primary statistics:</dt>
<dd>- <b>Hexagons:</b> the frequency of bins of position and copy number across the population being visualized.<br></dd>
<dd>- <b>Points:</b> the normalized-coverage for specific samples subset to only the intervals where the absolute value of the y-axis exceeds certain cutoffs.</dd><br>
<dt>This implementation removed intervals with no read dpeth, and thus omits homozygous mutations.</dt>  
</p>
</body>
</html>
')
index<-unlist(strsplit(index,"\n"))
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










