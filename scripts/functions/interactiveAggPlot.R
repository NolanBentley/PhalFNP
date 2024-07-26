aggPlotFun <- function(plottedDf,fileVec,genes,intMinDist = 1000000){
  library(ggplot2)
  library(hexbin)
  library(ggh4x)
  library(plotly)
  library(htmlwidgets)
  
  #Simplify
  plottedDf$IntervalMidpoint_Mbp <-round(plottedDf$IntervalMidpoint_Mbp,7)
  plottedDf$median <-round(plottedDf$median,3)
  genes$x   <-round(genes$x/1000000,6)
  genes$xend<-round(genes$xend/1000000,6)
  genes$y   <-round(genes$yOffset,6)
  genes$yend<-round(genes$yOffset,6)
  
  #Make a subset of interactive points
  fileVecSub   <- fileVec  [plottedDf$hasDivergentMedian ]
  plottedDfSub <- plottedDf[plottedDf$hasDivergentMedian,]
  
  #Loop across unique files
  uniFiles <- sort(unique(fileVecSub))
  
  #Load gene data
  genes$chr<-genes$seqid
  genes$yOffset <- genes$yOffset - min(genes$yOffset)
  genes$yOffset <- -genes$yOffset - 0.5
  
  i<-7
  for(i in 1:length(uniFiles)){
    currFile  <- uniFiles[i]
    currDf_i    <- plottedDf   [fileVec   ==currFile,]
    currDfSub_i <- plottedDfSub[fileVecSub==currFile,]
    uniChr    <- unique(plottedDf$chr[fileVec==currFile])
    
    #Modify if scaffold
    if(length(uniChr)&uniChr[1]=="scaffold"){
      currDf_i$chr      <- currDf_i$chr_orig
      currDfSub_i$chr   <- currDfSub_i$chr_orig
      uniChr            <- unique(plottedDf$chr_orig[fileVec==currFile])
    }
    
    #Modify for plotting
    currDf_i$id    <- as.factor(currDf_i$id)
    currDfSub_i$id <- factor(currDfSub_i$id,levels=levels(currDf_i$id))
    
    xLims <- c(min(c(0,currDf_i$IntervalMidpoint_Mbp)),max(currDf_i$IntervalMidpoint_Mbp))
    yLims <- c(min(c(0,genes$yOffset,currDf_i$median)),max(c(currDf_i$median)+0.1*(max(currDf_i$median)-min(currDf_i$median)),8))
    
    j<-1
    for(j in 1:length(uniChr)){
      currChr  <- uniChr[j]
      currDf   <- currDf_i   [currDf_i$chr==currChr   ,]
      currDfSub<- currDfSub_i[currDfSub_i$chr==currChr,]
      
      #Initialize plot
      currPlot <- ggplot(currDf,aes(IntervalMidpoint_Mbp,median,
                                    fill = round(..count../length(unique(currDf$id)),2),
                                    text = sprintf("Count per id: %0.5f", round(..count../length(unique(currDf$id)),2))))+
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
        #Subset to genes near divergent RD regions
        currGenes <- genes[genes$seqid%in%currChr,]
        bufferEdges <- intMinDist
        divPositions <- currDf[currDf$hasDivergentMedian,]
        divPositions <- divPositions[order(divPositions$chr_orig,divPositions$min,divPositions$max),]
        divPositions$breakPoints <- cumsum(c(1,diff(divPositions$min)>=intMinDist))
        divGMin <- aggregate(divPositions$min,by=list(divPositions$breakPoints),min)
        divGMax <- aggregate(divPositions$min,by=list(divPositions$breakPoints),max)
        currGenes$inDivG <- FALSE
        for(k in 1:nrow(divGMin)){
          currGenes$checkBeforeK <-(currGenes$start+bufferEdges)>(divGMin$x[k])
          currGenes$checkAfterK  <-(currGenes$end  -bufferEdges)<(divGMax$x[k])
          currGenes$inDivG <- currGenes$inDivG|(currGenes$checkBeforeK&currGenes$checkAfterK)
        }
        currGenes<- currGenes[currGenes$inDivG,]
        #Draw line segments
        if(nrow(currGenes)>0){
          currPlot <- currPlot+geom_segment(
            data = currGenes,
            mapping = aes(x=x,xend=xend,y=yOffset,yend=yOffset,fill=NULL,text=NULL),
            color="black"
          )
        }
        #Build hovertext
        currGenes$hoverText <- paste0(
          "<b>id:  </b>",currGenes$id," (",currGenes$strand," strand)","<br>",
          "<b>chr_range: </b>",trimws(currGenes$seqid),"_",
          gsub(" ","",paste0(format(currGenes$start,big.mark = ","),"-",format(currGenes$end,big.mark = ","))),
          "<br><b>other attributes:</br></b>",gsub(";","<br>",trimws(currGenes$attributes))
        )
        
        #Add label
        currPlot <- currPlot + 
          annotate("text",label="Genes:",color="black",
                   x=min(currGenes$x)-1.3,y=mean(currGenes$yOffset))
      }
      
      #Annotate name
      anno <- paste0(currChr,": ",gsub("LSVPlot_|\\.html","",basename(currFile)))
      currPlot <- currPlot + 
        annotate("text",label=anno,color="black",
                 x=mean(xLims),y=yLims[2]-0.1)
      
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
        #Muti chr single line plot
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
            x = c(currDfSub$IntervalMidpoint_Mbp,currGenes$x),
            y = c(currDfSub$median,currGenes$yOffset),
            color=as.factor(c(as.numeric(currDfSub$id)+1,rep(1,nrow(currGenes)))),
            fill="none",
            customdata = c(gsub("^\\.","/PhalFNP",currDfSub$singleLineMultiChr),currGenes$link),
            showlegend = F,
            hovertemplate= c(template,currGenes$hoverText)
          )
      }else{
        #Single chr single line plot
        currPlotly <- currPlotly %>%
          add_trace(
            type="scatter",mode = "markers",
            x = c(currDfSub$IntervalMidpoint_Mbp,currGenes$x),
            y = c(currDfSub$median,currGenes$yOffset),
            color=as.factor(c(as.numeric(currDfSub$id)+1,rep(1,nrow(currGenes)))),
            fill="none",
            customdata = c(gsub("^\\.","/PhalFNP",currDfSub$multiLineSingleChr),currGenes$link),
            showlegend = F,
            hovertemplate= c(template,currGenes$hoverText)
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
    saveWidget(currPlotly,file = currFile,selfcontained = F,title = gsub("LSVPlot_5x5Kbx40Kb_","",basename(currFile)))
    print(paste0(currFile," (",i," of ", length(uniFiles),")"))
  }
}