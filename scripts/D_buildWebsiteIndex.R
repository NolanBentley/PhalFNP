  ###Setup####
  #Set variables
  wd <- "~/Experiments/PhalFNP/"; setwd(wd)
  singlePrefix <- "phal_FIL.._"
  
  #Load functions
  source("~/Experiments/PhalFNP/scripts/functions/prepareDepths.R")
  
  #### Find files
  multiChrMultiLine   <- list.files("depthImages/multiChr/" ,full.names = T,pattern = "\\.html$")
  singleChrMultiLine  <- list.files("depthImages/singleChr/",full.names = T,pattern = "\\.html$")
  multiChrSingleLine  <- list.files("depthImages/multiChr/" ,full.names = T,pattern = paste0(singlePrefix,".*\\.html$"),recursive = T)
  singleChrSingleLine <- list.files("depthImages/singleChr/",full.names = T,pattern = paste0(singlePrefix,".*\\.html$"),recursive = T)
  
  #Make id bins
  singleChrSingleLine_bin <- paste0(dirname(singleChrSingleLine),".html")
  multiChrSingleLine_bin  <- paste0(dirname(multiChrSingleLine),".html")
  
  
  
  
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
  This website houses interactive plots for exploring large (>0.1 Mb) copy-number variants from the <i>Panicum. hallii</i> fast neutron population (FNP).<br>
  <br>
  <dt>The y-axis values were calculated to approximate ploidy across the x-axis intervals:</dt>
  <dd>- Each observation reports copy number across 7 intervals of ~10Kb.</dt>
  <dd>- The x-axis position is the midpoint of the 7 intervals being summarized.</dt>
  <dd>- The y-axis is the median value of copy number across the 7 intervals.</dt>
  <br>
  <dt>This data is visualized with two primary statistics:</dt>
  <dd>- <b>Hexagons:</b> the frequency of bins of position and copy number across the population being visualized.<br></dd>
  <dd>- <b>Points:</b> the median value of copy number across intervals where this value exceeds 0.5 away from diploid.</dd>
  </p>
  </body>
  </html>
  ')
  index<-unlist(strsplit(index,"\n"))
  index_endBody <- which(trimws(index)=="</body>")
  singleLineBinPages <-paste0("./depthImages/singleLineIndexes/",gsub("\\.html$","",basename(multiChrSingleLine_bin)),".html")
  index_added <- c(
    "<br><br><b><u>multi-sample multi-chromosome plots</b></u>",
    "<br>Clicking on points in these plots will take you to that chromosome's multi-sample plot.",
    "<br>Note: These are large (~40Mb) files and will take a bit to load and render.",
    "<ul>",
    sort(unique(multiChrMultiLine)),
    "</ul><br><br><b><u>multi-sample single-chromosome plots</b></u>",
    "<br>Clicking on points in these plots will take you to that sample's multi-chromosome plot. Clicking on those points will take you to the single-chromosome plot for that line.",
    "<ul>",
    sort(unique(singleChrMultiLine)),
    "</ul><br><br><b><u>Single line plots</b></u>",
    "This section takes you to single-sample plots based on the sample ID.",
    "<ul>",
    sort(unique(singleLineBinPages)),"</ul>"
  )
  
  #Add in links
  htmlPos  <-grep("html",index_added)
  links    <- gsub("^depthImages","/PhalFNP/depthImages",index_added[htmlPos])
  linkText <- gsub("^.*_|\\.html$","",basename(index_added[htmlPos]))
  linkText[linkText=="Chr"]<-"All chromosomes (This might cause your browser to become non-responsive!)"
  index_mod <- index_added
  index_mod[htmlPos] <- paste0('<li><a href="',links,'">',linkText,"</a></li>")
  
  #Put it all together
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
    currBin <- gsub("\\.html","",basename(posHtmlFiles)[i])
    index <- c(
      '<!DOCTYPE html>','<html>','<head>','<title>',currBin,'</title>','</head>',
      '<body>','<h1>','Welcome to the large variant explorer for ',currBin,'</h1>',
      '</body>','</html>'
    )
    index_endBody <- which(trimws(index)=="</body>")
    index_added <- c(
      "<br><br><b><u>single-sample multi-chromosome plots</b></u>",
      "<br>Clicking on non-diploid points in these plots will take you to that chromosome's single-sample plot.",
      "<ul>",
      sort(unique(multiChrSingleLine[basename(dirname(multiChrSingleLine))%in%currBin])),
      "</ul><br><br><b><u>single-sample single-chromosome plots</b></u>",
      "<br>Clicking on points in these plots will take you to that sample's multi-chromosome plot. Clicking on those points will take you to the single-chromosome plot for that line.",
      "<ul>",
      sort(unique(singleChrSingleLine[basename(dirname(singleChrSingleLine))%in%currBin]))
    )
    htmlPos  <-grep("html",index_added)
    links    <- gsub("^depthImages","/PhalFNP/depthImages",index_added[htmlPos])
    linkText <- gsub("LSVPlot_7x10kb_|phal_|\\.html$","",basename(index_added[htmlPos]))
    index_mod <- index_added
    index_mod[htmlPos] <- paste0('<li><a href="',links,'">',linkText,"</a></li>")
    
    #Add breaks between FIL30 and FIl20 and Chr# 
    FilNum <- gsub(".*FIL(..).*","\\1",index_mod[htmlPos])
    ChrNum <- gsub(".*Chr(..).*","\\1",index_mod[htmlPos])
    addLine <- c(abs(diff(as.numeric(as.factor(FilNum))))+abs(diff(as.numeric(as.factor(ChrNum)))),0)>0
    index_mod[htmlPos][addLine] <- paste0(index_mod[htmlPos][addLine],"<br>")
    
    index_built <- c(
      index[1:(index_endBody-1)],
      index_mod,
      index[index_endBody:length(index)]
    )
    write(index_built,posHtmlFiles[i])
  }
  
  
  
  
  
  
  
  
  
  
