
aggPlotFun <- function(x){
  
  #Create gradient
  rgbToHex <- function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}; apply(colorRamp(c("grey50","grey80"),space = "rgb")(0.1),1,rgbToHex)
  colsGrad <- apply(colorRamp(c("grey90","grey30"),space = "rgb")(log10(seq(1,10^1,length.out=100))),1,rgbToHex)
  
  #Initialize plot
  currPlot <- ggplot(currDf,aes(xColumn,yColumn,color=colorColumn))
  currPlot <- currPlot +
    geom_hex(data = currDf_orig,binwidth=c(1,1),mapping=aes(color=NULL))+
    scale_fill_gradientn(colours = colsGrad)
  currPlot <- currPlot + #geom_point(shape=21)+
    theme_bw()+
    guides(color="none")+
    facet_wrap(~chr_orig,nrow = 3)+
    labs(x=xLabel,y=yLabel,title=basename(currFile))
  currPlotly <- ggplotly(currPlot)
    currPlotly <- currPlotly %>% add_trace(
      type="scatter",mode = "markers",
      x = currDf[[xColumn]],
      y = currDf[[yColumn]],
      color=currDf[[colorColumn]],
      customdata = currDf$singlLineSinglChr,
      hovertemplate=paste(
        "<b>id:  </b>",currDf$id,"<br>",
        "<b>chr: </b>",currDf$chr_orig,"<br>",
        "<b>Interval (Mbp): </b>",currDf$min/1000000,"-",currDf$max/1000000,"<br>",
        "<b>Median value: </b>",currDf$median,"<br>"
      )
    )
    
    
}