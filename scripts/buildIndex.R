#### Load aggDf ####
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
  
  '<br><br>The y-axis is calculated from the median value across 5 intervals. ',
  'Each interval is summarized by the normalized mean read depth across a 2.5-5kb range of positions excluding positions with 0 depth. ',
  '<br>Each mean read depth has been tranformed by: ',
  '<b>first)</b> subtracting the median mean read depth across that sample. ',
  '<b>second)</b> dividing by the sample\'s mad.',
  '<b>third)</b> substracting by the loci\'s median normalized value.',
  '</p>','</body>','</html>'
)
index_endBody <- which(trimws(index)=="</body>")
singleLineBinPages <- paste0("./depthImages/singleLineIndexes/",sort(unique(aggDf$idNum_bin)),".html")
index_added <- c(
  "<br><br><b><u>multi-sample multi-chromosome plots</b></u>",
  "<br>Clicking on points in these plots will take you to that chromosome's multi-sample plot.",
  "<br>Note: These are large (~40Mb) files and will take a bit to load and render.",
  "<ul>",
  sort(unique(aggDf$multiLineMultiChr)),
  "</ul><br><br><b><u>multi-sample single-chromosome plots</b></u>",
  "<br>Clicking on points in these plots will take you to that sample's single-chromosome plot (see next section)",
  "<ul>",
  sort(unique(aggDf$multiLineSingleChr)),
  "</ul><br><br><b><u>Single line plots</b></u>",
  "This section is still under development. This section takes you to single-sample plots based on the sample ID.",
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
