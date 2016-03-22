library(shiny)
library("ggplot2")
library("plyr")
library("dplyr")
library("DT")

haplo.sum <- read.table("satrovirens02102016_haplo_filter.tbl", stringsAsFactors = FALSE) %>% 
  tbl_df
colnames(haplo.sum) <- c("id", "locus", "haplo", "depth", "logP.call", "logP.miscall")

haplo.cutoff <- haplo.sum %>%
  group_by(locus, id) %>%
  summarise(hapl.three.pl.st = ifelse(length(depth) > 2, 0, 0),
            hapl.three.pl.end = ifelse(length(depth) > 2, sort(depth, decr=T)[3]-1, 0),
            hapl.one.st = ifelse(sum(depth==max(depth))==1 && length(depth) > 1, sort(depth, decr=T)[2], 0),
            hapl.one.end = ifelse(sum(depth==max(depth))==1, max(depth),0))


shinyServer(function(input, output) {
  
  ranges <- reactiveValues(y = NULL)
  rangesH <- reactiveValues(y = NULL)

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
      scale_x_log10()+
      theme_bw()+
      xlab("read coverage cutoff")+
      ylab("Individual ID")+
      scale_color_manual(name= "Haplotypes:", values=c("1"="light grey","2"= "#4BBA82", "3+"="#A48A82"))+
      theme(legend.position="bottom")+
      coord_cartesian(ylim=ranges$y)
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
      scale_size_continuous("Read Depth of the most common haplotype (log 10)")+
      theme_bw()+
      ylab("Individual ID")+
      xlab ("Depth Ratio of the second common haplotype : first common haplotype")+
      theme(legend.position="bottom")+
      xlim(c(0,1))+
      coord_cartesian(ylim=ranges$y)
  })
  
  
  output$haplDensityPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > input$coverageMin) 
      
    if (input$selectLocus != "ALL") 
      haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus) 
    
    if (input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 

    haplo.ct <- haplo.filter %>%
      group_by(locus, id) %>% 
      summarise(tot.hapl = n()) 
    
    haplo.tot.tbl <- haplo.ct %>% 
      group_by(locus, tot.hapl) %>%
      summarise(ct = n()) %>%
      ungroup() %>%
      group_by(locus) %>%
      mutate(frac = ct/sum(ct))
    
    ggplot()+ 
      geom_point(data=haplo.tot.tbl, aes(x = tot.hapl, y = locus, size=frac, color=frac))+
      #scale_x_log10()+
      xlab("total number of unique haplotypes in an individual")+
      ylab("Locus ID")+
      scale_color_continuous("fraction")+
      scale_size_continuous(guide=FALSE)+
      theme_bw()+
      theme(legend.position="bottom")+
      coord_cartesian(ylim=rangesH$y)
  })
  
  
  observeEvent(input$plot_dblclick, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      ranges$y <- c(brush$ymin, brush$ymax)      
    } else {
      ranges$y <- NULL
    }
  })
  
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$y <- c(brush$ymin, brush$ymax)      
    } else {
      ranges$y <- NULL
    }
  })
  
  observeEvent(input$plotH_dblclick, {
    brush <- input$plotH_brush
    if (!is.null(brush)) {
      rangesH$y <- c(brush$ymin, brush$ymax)      
    } else {
      rangesH$y <- NULL
    }
  })
  
  
  output$distPlot2 <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()   
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > input$coverageMin, locus == input$selectLocus)
    
    if (input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    haplo.profile.frac <- haplo.filter %>%
      group_by(haplo) %>%
      summarise(n=n()) %>%
      ungroup() %>%
      mutate(frac = n/sum(n)) 
    
    haplo.split.profile <- sapply(1:nrow(haplo.profile.frac), function(i) {
      char.split <- strsplit(haplo.profile.frac[i,]$haplo, "")
      n.char <- length(char.split[[1]])
      sapply(1:n.char, function(j) c(i, j, char.split[[1]][j], haplo.profile.frac[i,]$frac)) 
    }) %>% 
      matrix(., ncol=4, byrow=T) %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      tbl_df()
    
    colnames(haplo.split.profile) <- c("group", "pos", "seq", "frac")
    haplo.split.profile <- haplo.split.profile %>% mutate(pos=as.numeric(pos),
                                                          frac=as.numeric(frac),
                                                          group=as.numeric(group))
    
    
    g <- ggplot(haplo.split.profile, aes(x=factor(pos), y=seq, group=group, size=frac, color=factor(group))) +
      xlab("relative position")+
      ylab("sequence")
    
    if (length(unique(haplo.split.profile$pos))==1) {
      g <- g + geom_point(alpha=0.9)
    }
    else {
      g <- g + geom_path(alpha=0.9)
    }
    
    
      g + scale_size_continuous(guide=FALSE)+
      theme_bw()+
      theme(legend.position="bottom")  
  })
  
  output$distPlot3 <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()   
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > input$coverageMin, locus == input$selectLocus)
    
    if (input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    haplo.filter <- haplo.filter %>% 
      group_by(locus, id) %>%
      summarise(tot.depth = sum(depth))
    
    ggplot(haplo.filter, aes(x=locus, y=tot.depth)) +
      xlab("locus id")+
      ylab("total read depth per indiv")+
      geom_violin()+
      theme_bw()+
      scale_y_log10()
    
  })
  
  output$distPlot4 <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()   
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > input$coverageMin, locus == input$selectLocus)
    
    if (input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    ggplot(haplo.filter, aes(x=factor(id), y=depth)) +
      xlab("indiv id")+
      ylab("Read depth")+
      theme_bw()+
      geom_violin()+
      scale_y_log10()
    
  })
    
  output$haploTbl <- DT::renderDataTable({
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > input$coverageMin) %>%
      select(id, locus, haplo, depth)
    
    if (!is.null(input$selectLocus) && input$selectLocus != "ALL") 
      haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus) 
    
    if (!is.null(input$selectIndiv) && input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    DT::datatable(
      haplo.filter, options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )
    )
  })
  
  
})
