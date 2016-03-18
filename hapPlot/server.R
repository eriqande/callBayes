library(shiny)
library("ggplot2")
library("plyr")
library("dplyr")

haplo.sum <- read.table("../data/hap/tot.sum", stringsAsFactors = FALSE) %>% 
  tbl_df
colnames(haplo.sum) <- c("id", "locus", "haplo", "depth", "logP.call", "logP.miscall")



# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  output$uiLocus <- renderUI({
    selectInput("selectLocus", label = "Locus:", 
                sort(c(unique(haplo.sum$locus), "ALL")),
                selected = "ALL"
    )
  })
    output$uiIndiv <- renderUI({
      selectInput("selectIndiv", label = "Indiv:", 
                  sort(c(unique(haplo.sum$id), "ALL")),
                  selected = "ALL"
      )
   # if (is.null(input$input_type))
  #    return()
    
  })
  
  output$distPlot <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()
    
    haplo.sample <- haplo.cutoff %>% filter(locus== input$selectLocus)
    if (input$selectIndiv != "ALL")
      haplo.sample <- haplo.cutoff %>% filter(locus== input$selectLocus, id == input$selectIndiv)
      if (dim(haplo.sample)[1]==0)
        return ()
      
    
    ggplot()+ 
      geom_segment(data=haplo.sample, aes(x = hapl.one.st, xend = hapl.one.end, y = id, yend = id, colour= "1"), size=2 )+
      geom_segment(data=haplo.sample, aes(x = hapl.three.pl.end, xend = hapl.one.st, y = id, yend = id, colour="2"), size=2 )+
      geom_segment(data=haplo.sample, aes(x = hapl.three.pl.st, xend = hapl.three.pl.end, y = id, yend = id, colour="3+"), size=1)+
      #scale_x_log10()+
      theme_bw()+
      xlab("read coverage cutoff")+
      ylab("Individual ID")+
      scale_color_manual(name= "Haplotypes:", values=c("1"="light grey","2"= "#4BBA82", "3+"="#A48A82"))+
      theme(legend.position="bottom")
  })
  
  output$distPlot1 <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()
    
    haplo.sample <- haplo.cutoff %>% filter(locus== input$selectLocus)
    if (input$selectIndiv != "ALL")
      haplo.sample <- haplo.cutoff %>% filter(locus== input$selectLocus, id == input$selectIndiv)
    if (dim(haplo.sample)[1]==0)
      return ()
    
    ggplot(data=haplo.sample, aes(x=hapl.one.st/hapl.one.end, y=id,size=log(hapl.one.end, 10)))+
      geom_point()+
      scale_size_continuous("Read Depth of the most common haplotype(log 10)")+
      theme_bw()+
      ylab("Individual ID")+
      xlab ("Depth Ratio of the second common haplotype : first common haplotype")+
      theme(legend.position="bottom")+
      xlim(c(0,1))
  })
})
