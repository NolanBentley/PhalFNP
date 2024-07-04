aggPlotFun <- function(plottedDf,fileVec){
  library(ggplot2)
  library(hexbin)
  library(plotly)
  library(htmlwidgets)
  
  #Make a subset of interactive points
  fileVecSub   <- fileVec  [plottedDf$hasDivergentMedian ]
  plottedDfSub <- plottedDf[plottedDf$hasDivergentMedian,]
  
  #Create gradient
  ## Made because of how plotly parses the data
  rgbToHex <- function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}
  colsGrad <- apply(colorRamp(c("grey90","grey30"),space = "rgb")(log10(seq(1,10^1,length.out=100))),1,rgbToHex)
  #cat(paste0("\nc('",paste0(colsGrad,collapse = "','"),"')\n\n")) #c('#E5E5E5','#DFDFDF','#D9D9D9','#D5D5D5','#D0D0D0','#CCCCCC','#C8C8C8','#C4C4C4','#C0C0C0','#BDBDBD','#BABABA','#B7B7B7','#B4B4B4','#B1B1B1','#AEAEAE','#ACACAC','#A9A9A9','#A7A7A7','#A5A5A5','#A2A2A2','#A0A0A0','#9E9E9E','#9C9C9C','#9A9A9A','#989898','#969696','#949494','#939393','#919191','#8F8F8F','#8E8E8E','#8C8C8C','#8B8B8B','#898989','#888888','#868686','#858585','#838383','#828282','#818181','#7F7F7F','#7E7E7E','#7D7D7D','#7B7B7B','#7A7A7A','#797979','#787878','#777777','#767676','#757575','#737373','#727272','#717171','#707070','#6F6F6F','#6E6E6E','#6D6D6D','#6C6C6C','#6B6B6B','#6A6A6A','#696969','#686868','#686868','#676767','#666666','#656565','#646464','#636363','#626262','#626262','#616161','#606060','#5F5F5F','#5E5E5E','#5E5E5E','#5D5D5D','#5C5C5C','#5B5B5B','#5A5A5A','#5A5A5A','#595959','#585858','#585858','#575757','#565656','#555555','#555555','#545454','#535353','#535353','#525252','#515151','#515151','#505050','#505050','#4F4F4F','#4E4E4E','#4E4E4E','#4D4D4D','#4D4D4D')
  
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
        labs(x="Interval midpoint (Mbp)",y="Median(average read depth)")+
        coord_cartesian(xlim=xLims,ylim=yLims)
      
      #Add hex plot
      currPlot <- currPlot +
        geom_hex(binwidth=c(1,1))+
        scale_fill_gradientn(colours = colsGrad)
      
      #Annotate name
      anno <- paste0(currChr,": ",gsub("LSVPlot_|\\.html","",basename(currFile)))
      currPlot <- currPlot + 
        annotate("text",label=anno,color="black",
                 x=mean(xLims),y=yLims[2]+0*(yLims[2]-yLims[1]))
      
      #currPlotly <- currPlotly %>% 
      # layout(annotations = list(x=0.1,y=1.03,text=anno,showarrow=F))#,xref="paper",yref="paper"))
      
      #Add interactive scatterplot
      currPlotly <- ggplotly(currPlot) %>%
        add_trace(
          type="scatter",mode = "markers",
          x = currDfSub$IntervalMidpoint_Mbp,
          y = currDfSub$median,
          color=currDfSub$id,
          customdata = gsub("^\\.","/PhalFNP",currDfSub$singleLineSingleChr),
          showlegend = F,
          hovertemplate=paste(
            "<b>id:  </b>",currDfSub$id,"<br>",
            "<b>chr: </b>",currDfSub$chr_orig,"<br>",
            "<b>Interval start (bp): </b>",format(currDfSub$min,big.mark = ","),"<br>",
            "<b>Interval end (bp): </b>",format(currDfSub$max,big.mark = ","),"<br>",
            "<b>Interval range (bp): </b>",format(currDfSub$max - currDfSub$min +1,big.mark = ","),"<br>",
            "<b>Median value: </b>",round(currDfSub$median,3)
          )
        )
        
  	  #Hyperlink script from: https://stackoverflow.com/questions/71819237/ggplotly-clickable-link-in-r-plot
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
      
      
      #Add to combined plots
      if(j==1){plotList <- list()}
      plotList[[j]] <- currPlotly
    }
    if(length(plotList)>3){
      currPlotly <- subplot(plotList,nrows = 3,shareX = T,shareY = T,titleX = T,titleY = T)
    }else{
      currPlotly <- subplot(plotList,nrows = length(plotList),shareX = T,shareY = T,titleX = T,titleY = T)
    }
    dir.create(dirname(currFile),recursive = T,showWarnings = F)
    saveWidget(currPlotly,file = currFile,selfcontained = T)
    print(paste0(currFile," (",i," of ", length(uniFiles),")"))
  }
}