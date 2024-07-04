aggPlotFun <- function(plottedDf,fileVec){
  library(ggplot2)
  library(hexbin)
  library(plotly)
  library(htmlwidgets)
  
  #Make a subset of interactive points
  fileVecSub   <- fileVec  [plottedDf$hasDivergentMedian ]
  plottedDfSub <- plottedDf[plottedDf$hasDivergentMedian,]
  
  #Loop across unique files
  uniFiles <- sort(unique(fileVecSub))
  for(i in 1:length(uniFiles)){
    currFile  <- uniFiles[i]
    currDf_i    <- plottedDf   [fileVec   ==currFile,]
    currDf_i$id <- as.factor(currDf_i$id)
    currDfSub_i <- plottedDfSub[fileVecSub==currFile,]
    currDfSub_i$id <- factor(currDfSub_i$id,levels=levels(currDf_i$id))
    
    uniChr    <- unique(plottedDf$chr[fileVec==currFile])
    xLims <- c(min(currDf_i$IntervalMidpoint_Mbp),max(currDf_i$IntervalMidpoint_Mbp))
    yLims <- c(min(currDf_i$median),max(currDf_i$median)+0.1*(max(currDf_i$median)-min(currDf_i$median)))
    for(j in 1:length(uniChr)){
      currChr  <-uniChr[j]
      currDf   <- currDf_i   [currDf_i$chr==currChr   ,]
      currDfSub<- currDfSub_i[currDfSub_i$chr==currChr,]
      
      #Initialize plot
      currPlot <- ggplot(currDf,aes(IntervalMidpoint_Mbp,median))+
        theme_bw()+
        guides(color="none")+
        labs(x="Interval midpoint (Mbp)\n ",y="\n Median(average read depth)")+
        coord_cartesian(xlim=xLims,ylim=yLims)
      
      #Add hex plot
      currPlot <- currPlot +
        geom_hex(binwidth=c(1,1))+
        scale_fill_gradient(low = "grey90",high = "grey30",trans="log10",name = "log10(count)")
      
      #Annotate name
      anno <- paste0(currChr,": ",gsub("LSVPlot_|\\.html","",basename(currFile)))
      currPlot <- currPlot + 
        annotate("text",label=anno,color="black",
                 x=mean(xLims),y=yLims[2]+0*(yLims[2]-yLims[1]))
      
      #Add interactive scatterplot
      currPlotly <- ggplotly(currPlot)
      
      template <- paste(
        "<b>id:  </b>",currDfSub$id,"<br>",
        "<b>chr: </b>",currDfSub$chr_orig,"<br>",
        "<b>Interval start (bp): </b>",format(currDfSub$min,big.mark = ","),"<br>",
        "<b>Interval end (bp): </b>",format(currDfSub$max,big.mark = ","),"<br>",
        "<b>Interval range (bp): </b>",format(currDfSub$max - currDfSub$min +1,big.mark = ","),"<br>",
        "<b>Median value: </b>",round(currDfSub$median,3)
      )
      if(length(uniChr)>1){
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = currDfSub$IntervalMidpoint_Mbp,
            y = currDfSub$median,
            color=currDfSub$id,
            customdata = gsub("^\\.","/PhalFNP",currDfSub$multiLineSingleChr),
            showlegend = F,
            hovertemplate= template
          )
      }else if(length(unique(currDfSub_i$id))>1){
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = currDfSub$IntervalMidpoint_Mbp,
            y = currDfSub$median,
            color=currDfSub$id,
            customdata = gsub("^\\.","/PhalFNP",currDfSub$singleLineMultiChr),
            showlegend = F,
            hovertemplate= template
          )
      }else{
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = currDfSub$IntervalMidpoint_Mbp,
            y = currDfSub$median,
            color=currDfSub$id,
            customdata = gsub("^\\.","/PhalFNP",currDfSub$multiLineSingleChr),
            showlegend = F,
            hovertemplate= template
          )
      }
        
      #Add to combined plots
      if(j==1){plotList <- list()}
      plotList[[j]] <- currPlotly
    }
    #Combine if needed and add hyperlinks
    if(length(plotList)>1){
      if(length(plotList)>2){
        currPlotly <- subplot(plotList,nrows = 3,shareX = T,shareY = T,titleX = T,titleY = T)
      }else{
        currPlotly <- subplot(plotList,nrows = length(plotList),shareX = T,shareY = T,titleX = T,titleY = T)
      }
    }
    currPlotly <- onRender(
        currPlotly, "
        function(el) {
          el.on('plotly_click', function(d) {
            var url = d.points[0].customdata;
            window.open(url);
          });
        }
        "
    )
    
    #Save output
    dir.create(dirname(currFile),recursive = T,showWarnings = F)
    saveWidget(currPlotly,file = currFile,selfcontained = T)
    print(paste0(currFile," (",i," of ", length(uniFiles),")"))
  }
}