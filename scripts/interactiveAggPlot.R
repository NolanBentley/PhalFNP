aggPlotFun <- function(plottedDf,fileVec,genes){
  library(ggplot2)
  library(hexbin)
  library(ggh4x)
  library(plotly)
  library(htmlwidgets)
  
  #Make a subset of interactive points
  fileVecSub   <- fileVec  [plottedDf$hasDivergentMedian ]
  plottedDfSub <- plottedDf[plottedDf$hasDivergentMedian,]
  
  #Loop across unique files
  uniFiles <- sort(unique(fileVecSub))
  
  #Load gene data
  genes$chr<-genes$seqid
  genes$yOffset <- genes$yOffset - min(genes$yOffset)
  
  for(i in 1:length(uniFiles)){
    currFile  <- uniFiles[i]
    currDf_i    <- plottedDf   [fileVec   ==currFile,]
    currDf_i$id <- as.factor(currDf_i$id)
    currDfSub_i <- plottedDfSub[fileVecSub==currFile,]
    currDfSub_i$id <- factor(currDfSub_i$id,levels=levels(currDf_i$id))
    
    uniChr    <- unique(plottedDf$chr[fileVec==currFile])
    xLims <- c(min(c(0,currDf_i$IntervalMidpoint_Mbp)),max(currDf_i$IntervalMidpoint_Mbp))
    yLims <- c(min(c(0,currDf_i$median)),max(c(currDf_i$median)+0.1*(max(currDf_i$median)-min(currDf_i$median)),8))
    for(j in 1:length(uniChr)){
      currChr  <-uniChr[j]
      currDf   <- currDf_i   [currDf_i$chr==currChr   ,]
      currDfSub<- currDfSub_i[currDfSub_i$chr==currChr,]
      
      #Initialize plot
      currPlot <- ggplot(currDf,aes(IntervalMidpoint_Mbp,median,
                                    fill = ..count../length(unique(currDf$id)),
                                    text = sprintf("Count per id: %0.5f", ..count../length(unique(currDf$id)))))+
        theme_bw()+
        guides(color="none")+
        labs(x="Interval midpoint (Mbp)\n ",y="\n Median(average read depth)")+
        scale_x_continuous(breaks     =seq(0,xLims[2],by=2),
                           minor_breaks=seq(0,xLims[2],by=0.02),
                           guide="axis_minor")+
        coord_cartesian(xlim=xLims,ylim=yLims)+
        theme(ggh4x.axis.ticks.length.minor=rel(1))
      
      #Add hex plot
      currPlot <- currPlot +
        geom_hex(binwidth=c(0.8,0.1))+
        scale_fill_gradient(low = "grey90",high = "grey50",name="Count\nper ID")
      
      #If only one chromosome, add in genes
      if(length(uniChr)==1){
        currPlot <- currPlot+geom_segment(
          data = genes[genes$seqid%in%currChr,],
          mapping = aes(x=x/1000000,xend=xend/1000000,y=yOffset,yend=yOffset,fill=NULL,text=attributes),
          color="black"
        )+ geom_point(
          data = genes[genes$seqid%in%currChr,],
          mapping = aes(x=x/1000000,y=yOffset,fill=NULL,text=attributes),
          arrow = arrow(length=unit(0.001,"cm")),color="black"
        )
      }
      
      #Annotate name
      anno <- paste0(currChr,": ",gsub("LSVPlot_|\\.html","",basename(currFile)))
      currPlot <- currPlot + 
        annotate("text",label=anno,color="black",
                 x=mean(xLims),y=yLims[2]+0*(yLims[2]-yLims[1]))
      
      #Add interactive scatterplot
      currPlotly <- ggplotly(currPlot,tooltip = "text")
      
      template <- paste(
        "<b>id:  </b>",currDfSub$id,"<br>",
        "<b>chr: </b>",currDfSub$chr_orig,"<br>",
        "<b>Interval start (bp): </b>",format(currDfSub$min,big.mark = ","),"<br>",
        "<b>Interval end (bp): </b>",format(currDfSub$max,big.mark = ","),"<br>",
        "<b>Interval range (bp): </b>",format(currDfSub$max - currDfSub$min +1,big.mark = ","),"<br>",
        "<b>Median value: </b>",round(currDfSub$median,3)
      )
      
      if(length(uniChr)>1&length(unique(currDfSub_i$id))>1){
        #Multi chr plot
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = currDfSub$IntervalMidpoint_Mbp,
            y = currDfSub$median,
            color=currDfSub$id,
            fill="none",
            customdata = gsub("^\\.","/PhalFNP",currDfSub$multiLineSingleChr),
            showlegend = F,
            hovertemplate= template
          )
      }else if(length(uniChr)>1&length(unique(currDfSub_i$id))==1){
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = currDfSub$IntervalMidpoint_Mbp,
            y = currDfSub$median,
            color=currDfSub$id,
            fill="none",
            customdata = gsub("^\\.","/PhalFNP",currDfSub$singleLineSingleChr),
            showlegend = F,
            hovertemplate= template
          )
      }else if(length(unique(currDfSub_i$id))>1){
        #Single chr multi line plot
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = currDfSub$IntervalMidpoint_Mbp,
            y = currDfSub$median,
            color=currDfSub$id,
            fill="none",
            customdata = gsub("^\\.","/PhalFNP",currDfSub$singleLineMultiChr),
            showlegend = F,
            hovertemplate= template
          )
      }else{
        #Single chr single line plot
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = currDfSub$IntervalMidpoint_Mbp,
            y = currDfSub$median,
            color=currDfSub$id,
            fill="none",
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
    saveWidget(currPlotly,file = currFile,selfcontained = F)
    print(paste0(currFile," (",i," of ", length(uniFiles),")"))
  }
}