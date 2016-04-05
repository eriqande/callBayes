library(shiny)
library("ggplot2")
library("plyr")
library("dplyr")
library("tidyr")
library("DT")
library("grid")
library("scales")

haplo.sum<- readRDS("satrovirens01092016_panel1_haplo_filter.rds")
#haplo.sum<- readRDS("satrovirens02102016_panel2_haplo_filter.rds")
#colnames(haplo.sum) <- c("id", "locus", "haplo", "depth", "logP.call", "logP.miscall", "allele.balance","rank")
haplo.sum <- haplo.sum %>% mutate(id = as.character(id))

# haplo.cutoff <- haplo.sum %>%
#   group_by(locus, id) %>%
#   summarise(hapl.three.pl.st = ifelse(length(depth) > 2, 0, 0),
#             hapl.three.pl.end = ifelse(length(depth) > 2, sort(depth, decr=T)[3]-1, 0),
#             hapl.one.st = ifelse(sum(depth==max(depth))==1 && length(depth) > 1, sort(depth, decr=T)[2], 0),
#             hapl.one.end = ifelse(sum(depth==max(depth))==1, max(depth),0))

n.locus <- length(unique(haplo.sum$locus))
n.indiv <- length(unique(haplo.sum$id))
locus.label.tbl <-  data.frame(locus =sort(unique(haplo.sum$locus)), stringsAsFactors = F) %>% tbl_df()
locus.label <- c("ALL",sort(unique(haplo.sum$locus)))
indiv.label.tbl <-  data.frame(id =sort(unique(haplo.sum$id)), stringsAsFactors = F) %>% tbl_df()
indiv.label <- c("ALL",sort(unique(haplo.sum$id)))

shinyServer(function(input, output, session) {
  
  ranges <- reactiveValues(y = NULL, x = NULL)
  rangesH <- reactiveValues(y = NULL)
  filterParam <- reactiveValues(minRead = 1, minAllele = 0.2)
  plotParam<- reactiveValues(byLocus.width= length(locus.label)*9)
  
  ## updating Locus and individidual choice at the start of the session:
  updateSelectInput(session, "selectLocus", selected="ALL", choices=locus.label)
  updateSelectInput(session, "selectIndiv", selected="ALL", choices=indiv.label)
  
  # reacting to the locus & Indiv's previous and next button 
  observeEvent(input$locusBack, {
    indx <- isolate(which(locus.label==input$selectLocus))
    label <- ifelse(indx > 1, locus.label[indx-1], locus.label[indx])
    updateSelectInput(session, "selectLocus", selected=label)
  })
  observeEvent(input$locusFor, {
    indx <- isolate(which(locus.label==input$selectLocus))
    label <- ifelse(indx < length(locus.label), locus.label[indx+1], locus.label[indx])     
    updateSelectInput(session, "selectLocus", selected=label)
  })
  observeEvent(input$indivBack, {
    indx <- isolate(which(indiv.label==input$selectIndiv))
    label <- ifelse(indx > 1, indiv.label[indx-1], indiv.label[indx])
    updateSelectInput(session, "selectIndiv", selected=label)
  })
  observeEvent(input$indivFor, {
    indx <- isolate(which(indiv.label==input$selectIndiv))
    label <- ifelse(indx < length(indiv.label), indiv.label[indx+1], indiv.label[indx])     
    updateSelectInput(session, "selectIndiv", selected=label)
  })
  
  # reacting to the filter update button
  observeEvent(input$updateFilter, {
    filterParam$minRead <- input$coverageMin
    filterParam$minAllele <- input$minAlleleRatio
  })
  
  
  haplo.summaryTbl <- reactive({
    haplo.filter <- haplo.sum %>% 
      filter(depth > filterParam$minRead, rank <= 2, allele.balance >= filterParam$minAllele) 
    
    if (input$selectLocus != "ALL") haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus)
    if (input$selectIndiv != "ALL") haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    haplo.filter <- haplo.filter %>% 
      group_by(locus, id) %>%
      arrange(-depth) %>%
      summarise(haplotype.1 = ifelse(length(depth)==1, haplo[1], sort(haplo)[1]),
                haplotype.2 = ifelse(length(depth)==1, haplo[1], sort(haplo)[2]),
                read.depth.1 = depth[1],
                read.depth.2 = ifelse(length(depth)==1, depth[1], depth[2])) 
  })
  
  haplo.freqTbl <- reactive({
    
    obs.freq.tbl<-  haplo.summaryTbl() %>%
      ungroup() %>%
      group_by(locus) %>% 
      mutate(tot.haplo = n()) %>%
      group_by(locus, haplotype.1, haplotype.2) %>%
      summarise(obs.freq=n()/tot.haplo[1])
    
    expect.freq.tbl <- gather(obs.freq.tbl, whichHap,  hap1, 2:3) %>%
      group_by(locus, hap1) %>% 
      summarise(n=sum(obs.freq/2)) %>%
      mutate(hap2 = hap1, n1 = n) %>%
      expand(., nesting(hap1, n), nesting(hap2, n1)) %>%
      mutate(expected.freq = ifelse(hap1==hap2, n*n1, 2*n*n1)) %>%
      group_by(locus, hap1, hap2) %>%
      mutate(haplotype.1 = sort(c(hap1, hap2))[1],
             haplotype.2 = sort(c(hap1, hap2))[2]) %>%
      ungroup() %>%
      select(locus, haplotype.1, haplotype.2, expected.freq) %>%
      distinct()
      
    inner_join(obs.freq.tbl, expect.freq.tbl, by=c("locus", "haplotype.1", "haplotype.2"))
  })
  
  output$downloadData <- downloadHandler(
    filename = 'filtered_haplotype.csv',
    content = function(file) {    
      write.csv(haplo.summaryTbl() %>% 
                  rename("Indiv.ID"=id),
                file)
    }
  )
  
  Filter.haplo.sum <- reactive({
    haplo.filter <- haplo.sum %>% 
      filter(depth > filterParam$minRead, allele.balance >= filterParam$minAllele) 
    
    if(input$topTwo)
      haplo.filter <- haplo.filter %>% filter(rank <= 2)
    
    if (input$selectLocus != "ALL") 
      haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus) 
    
    if (input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    haplo.filter
  })
  
  Get.tbl.by.locus <- reactive({
    haplo.ct <- Filter.haplo.sum() %>%
      group_by(locus, id) %>% 
      summarise(tot.hapl = n(), tot.depth = sum(depth))     
  })
  
  Get.tbl.by.id <- reactive({
    haplo.ct <- Filter.haplo.sum() %>%
      group_by(id, locus) %>% 
      summarise(tot.depth = sum(depth))     
  })
  
  
  
  # BY LOCUS PANEL::  
  
  output$haplDensityPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    haplo.tot.tbl <- Get.tbl.by.locus() %>% 
      group_by(locus, tot.hapl) %>%
      summarise(ct = n()) %>%
      ungroup() %>%
      group_by(locus) %>%
      mutate(frac = ct/sum(ct))
    
    uniqH.perI.tbl <- right_join(haplo.tot.tbl, locus.label.tbl, by="locus") 
    uniqH.perI.tbl[is.na(uniqH.perI.tbl)]<- 0
    
    if (input$selectLocus != "ALL") {
      uniqH.perI.tbl <- uniqH.perI.tbl %>% filter(locus == input$selectLocus)
    }
    
    ggplot()+ 
      geom_point(data=uniqH.perI.tbl, aes(x = tot.hapl, y = locus, size=frac, color=frac))+
      #scale_x_log10()+
      xlab("number of unique haplotypes per individual")+
      ylab("Locus ID")+
      scale_color_continuous(guide=FALSE)+#"fraction")+
      scale_size_continuous(guide=FALSE)+
      theme_bw()+
      theme(legend.position="bottom",
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      #scale_x_discrete(breaks= pretty_breaks())+
      coord_cartesian(ylim=rangesH$y)
  },height = function(){max(ifelse(input$selectLocus=="ALL",9*length(locus.label),1),400) })  
  
  output$numHapPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    
    frac.calleable <- haplo.summaryTbl() %>% group_by(locus) %>% summarise(n=length(unique(c(haplotype.1,haplotype.2))))
    
    frac.calleable <- right_join(frac.calleable, locus.label.tbl, by="locus") 
    frac.calleable[is.na(frac.calleable)]<- 0
    
    if (input$selectLocus != "ALL") {
      frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    }
    
    
    ggplot(frac.calleable, aes(x=n, y=locus))+
      geom_point()+
      xlab("num of overall unique haplotypes")+
      ylab("")+
      scale_size_continuous(guide=FALSE)+#"fraction")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 2, 0, 0), "mm"))+
      coord_cartesian(ylim=rangesH$y)#+
      #scale_x_discrete(limits=c(-1, max(frac.calleable$n)+1)) #breaks= pretty_breaks()
  },height = function(){max(ifelse(input$selectLocus=="ALL",9*length(locus.label),1),400) })
  
  
  output$fracIndivPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    frac.calleable <- haplo.summaryTbl() %>% group_by(locus) %>% summarise(f=n()/n.indiv)
    frac.calleable <- right_join(frac.calleable, locus.label.tbl, by="locus") 
    frac.calleable[is.na(frac.calleable)]<- 0
    
    if (input$selectLocus != "ALL") {
      frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    }
    
    ggplot(frac.calleable, aes(x=f, y=locus, color=f))+
      geom_point()+
      xlab("fraction of indiv w/ calleable haplotype")+
      ylab("")+
      scale_color_continuous(guide=FALSE)+#"fraction")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      coord_cartesian(ylim=rangesH$y)+
      xlim(c(0,1))
  },height = function(){max(ifelse(input$selectLocus=="ALL",9*length(locus.label),1),400) })
  
  output$readDepthPerLocus <- renderPlot({
    
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()   
    
    readDepth.perI.tbl <- right_join(Get.tbl.by.locus(), locus.label.tbl, by="locus") 
    readDepth.perI.tbl[is.na(readDepth.perI.tbl)]<- 0
    
    if (input$selectLocus != "ALL") {
      readDepth.perI.tbl <- readDepth.perI.tbl %>% filter(locus == input$selectLocus)
    }
    
    ggplot(readDepth.perI.tbl, aes(x=locus, y=tot.depth)) +
      xlab("")+
      ylab("total read depth per individual")+
      geom_violin()+
      #geom_point()+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      scale_y_log10()+
      coord_flip(xlim=rangesH$y)    
  },height = function(){max(ifelse(input$selectLocus=="ALL",9*length(locus.label),1),400) })
  
  ## BY INDIVIDUAL PANEL::   
  
  output$AlleleRatioByIndiv <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > filterParam$minRead, allele.balance >= filterParam$minAllele, rank <= 2) 
    
    if (input$selectLocus != "ALL") haplo.filter <- haplo.filter %>% filter(locus== input$selectLocus)   
    if (input$selectIndiv != "ALL") haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    
    if (dim(haplo.filter)[1]==0) return ()
    
    haplo.filter <- haplo.filter %>% 
      group_by(locus, id) %>% 
      summarise(depth.ratio = ifelse(length(depth)==1, 0, min(allele.balance)),
                depth.first = max(depth))
    
    haplo.filter <- right_join(haplo.filter, indiv.label.tbl, by="id") 
    haplo.filter[is.na(haplo.filter)]<- 0
    
    if (input$selectIndiv != "ALL") {
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    }
    
    ggplot(data=haplo.filter, aes(x=depth.ratio, y=id,size=log(depth.first, 10)))+
      geom_point(alpha=0.4)+
      scale_size_continuous(guide=FALSE)+#"Read Depth of the most common haplotype (log 10)")+
      theme_bw()+
      ylab("individual ID")+
      xlab ("depth ratio of the second common haplotype : first common haplotype")+
      theme(legend.position="bottom")+
      xlim(c(0,1))+
      coord_cartesian(ylim=ranges$y)+
      geom_vline(xintercept=filterParam$minAllele, linetype="dashed", color = "red")
  },height = function(){max(ifelse(input$selectIndiv=="ALL",9*length(indiv.label),1),400) })
  
  output$fracHaploPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    haplo.filter <- haplo.summaryTbl() %>% group_by(id) %>% summarise(f=n()/n.locus)
    haplo.filter <- right_join(haplo.filter, indiv.label.tbl, by="id") 
    haplo.filter[is.na(haplo.filter)]<- 0
    
    if (input$selectIndiv != "ALL") {
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    }
    
    ggplot(haplo.filter, aes(x=f, y=id, color=f))+
      geom_point()+
      xlab("fraction of calleable haplotypes")+
      ylab("")+
      scale_color_continuous(guide=FALSE)+#"fraction")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'))+
            #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      coord_cartesian(ylim=ranges$y)+
      xlim(c(0,1))
  },height = function(){max(ifelse(input$selectIndiv=="ALL",9*length(indiv.label),1),400) })  
  
  
  output$meanReadDepthByIndiv <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    haplo.filter <- Get.tbl.by.id() %>%
      ungroup()%>%
      group_by(id) %>%
      summarise(mean.depth = mean(tot.depth))
    
    haplo.filter <- right_join(haplo.filter, indiv.label.tbl, by="id") 
    haplo.filter[is.na(haplo.filter)]<- 0.0001
    
    if (input$selectIndiv != "ALL") {
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    }
    
    ggplot(haplo.filter, aes(y=id, x=mean.depth)) +
      geom_point()+
      ylab("")+
      xlab("mean locus read depth ")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'))+
      #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      scale_x_log10()+
      coord_cartesian(ylim=ranges$y)
    
  },height = function(){max(ifelse(input$selectIndiv=="ALL",9*length(indiv.label),1),400) })
  
  
  output$readDepthByIndiv <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()
    
    haplo.filter <- right_join( Filter.haplo.sum(), indiv.label.tbl, by="id") 
    haplo.filter[is.na(haplo.filter)]<- 0
    
    if (input$selectIndiv != "ALL") {
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    }
  
    ggplot(haplo.filter, aes(x=id, y=depth, group=id)) +
      xlab("")+
      ylab("haplotype read depth")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'))+
            #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      geom_violin()+
      scale_y_log10()+
      coord_flip(xlim=ranges$y)
    
  },height = function(){max(ifelse(input$selectIndiv=="ALL",9*length(indiv.label),1),400) })
  
#   output$distPlot <- renderPlot({
#     if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
#       return ()
# 
#     haplo.sample <- haplo.cutoff %>% filter(locus== input$selectLocus)
#     
#     if (input$selectIndiv != "ALL")
#       haplo.sample <- haplo.sample %>% filter(locus== input$selectLocus, id == input$selectIndiv)
#     if (dim(haplo.sample)[1]==0)
#       return ()
#     
#     ggplot()+ 
#       geom_segment(data=haplo.sample, aes(x = hapl.one.st, xend = hapl.one.end, y = id, yend = id, colour= "1"), size=2 )+
#       geom_segment(data=haplo.sample, aes(x = hapl.three.pl.end, xend = hapl.one.st, y = id, yend = id, colour="2"), size=2 )+
#       geom_segment(data=haplo.sample, aes(x = hapl.three.pl.st, xend = hapl.three.pl.end, y = id, yend = id, colour="3+"), size=1)+
#       scale_x_log10()+
#       theme_bw()+
#       xlab("read coverage cutoff")+
#       ylab("Individual ID")+
#       scale_color_manual(name= "Haplotypes:", values=c("1"="light grey","2"= "#4BBA82", "3+"="#A48A82"))+
#       theme(legend.position="bottom")+
#       coord_cartesian(ylim=ranges$y)
#   })
#   
  
  
  
  
  
  
  observeEvent(input$plot_dblclick, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      ranges$y <- c(brush$ymin, brush$ymax)
      ranges$x <- c(brush$xmin, brush$xmax)
    } else {
      ranges$y <- NULL
      ranges$x <- NULL
    }
  })
  
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$y <- c(brush$ymin, brush$ymax)
      ranges$x <- c(brush$xmin, brush$xmax)
    } else {
      ranges$y <- NULL
      ranges$x <- NULL
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
  
  
  
#ABOUT HAPLOTYPE distribution panel  
  output$hapSeq <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()   
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > filterParam$minRead, locus == input$selectLocus, allele.balance >= filterParam$minAllele) 
    
    if(input$topTwo)
      haplo.filter <- haplo.filter %>% filter(rank <= 2)
    
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
      xlab("variant position")+
      ylab("sequence")+
      scale_size_continuous(range = c(3,20), guide=FALSE)
    
    if (length(unique(haplo.split.profile$pos))==1) {
      g <- g + geom_point(alpha=0.9)
    }
    else {
      g <- g + geom_path(alpha=0.9)
    }
    
    
    g + scale_size_continuous(guide=FALSE)+
      scale_color_discrete(guide=FALSE)+
      theme_bw()+
      theme(legend.position="bottom")  
  },height = function(){ifelse(input$selectLocus=="ALL",0,400) })
  
  output$histHap <- renderPlot({
    
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()   
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > filterParam$minRead, locus == input$selectLocus, allele.balance >= filterParam$minAllele) 
    
    if(input$topTwo)
      haplo.filter <- haplo.filter %>% filter(rank <= 2)
    
    if (input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    haplo.filter <- haplo.filter %>% group_by(haplo) %>% summarise(f=n()/n.indiv)
    ggplot(haplo.filter, aes(y=haplo, x = f, color=factor(haplo))) +
      geom_point(size=4)+
      xlab("Fraction of Individuals")+
      ylab("haplotype")+
      theme_bw()+
      scale_color_discrete(guide=FALSE)
    
  },height = function(){ifelse(input$selectLocus=="ALL",0,400) })
  
  
  
  output$PairWiseHap <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv))
      return ()   
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > filterParam$minRead, locus == input$selectLocus, allele.balance >= filterParam$minAllele) 
    
    #if(input$topTwo) 
    haplo.filter <- haplo.filter %>% filter(rank <= 2)
    if (input$selectIndiv != "ALL") haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    haplo.filter <- haplo.filter %>% 
      group_by(locus, id) %>%
      arrange(-depth) %>%
      summarise(hap1 = ifelse(length(depth)==1, haplo[1], sort(haplo)[1]),
                hap2 = ifelse(length(depth)==1, haplo[1], sort(haplo)[2])) %>%
      ungroup() %>%
      group_by(locus, hap1, hap2) %>%
      summarise(n=n())
    
    n.hap <- 2*sum(haplo.filter$n)
    freq.hap <- gather(haplo.filter, whichHap,  hap, 2:3) %>%
      group_by(locus, hap) %>% 
      summarise(n=sum(n)/n.hap) %>%
      mutate(hap1 = hap, n1 = n) %>%
      expand(., nesting(hap, n), nesting(hap1, n1)) %>%
      mutate(freq = ifelse(hap==hap1, n*n1, 2*n*n1)) %>%
      rename("hap1"=hap, "hap2"=hap1) %>%
      group_by(locus, hap1, hap2) %>%
      mutate(re.hap1 = sort(c(hap1, hap2))[2],
             re.hap2 = sort(c(hap1, hap2))[1])
    
    ggplot(haplo.filter, aes(x=hap1, y=hap2, size=n, color=hap1==hap2))+
      geom_point()+
      xlab("haplotype 1")+
      ylab("haplotype 2")+
      geom_point(data=freq.hap, aes(x=re.hap1,y=re.hap2, size=freq*n.hap/2), shape=21, fill=NA, color="black")+
      scale_color_discrete(guide=FALSE)+
      scale_size_continuous(range = c(3,20), guide=FALSE)+
      theme_bw()
  },height = function(){ifelse(input$selectLocus=="ALL",0,400) })
  
  
  output$haploTbl <- DT::renderDataTable({
    
    haplo.filter <- haplo.sum %>% 
      filter(depth > filterParam$minRead) %>%
      select(id, locus, haplo, depth)
    
    if (!is.null(input$selectLocus) && input$selectLocus != "ALL") 
      haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus) 
    
    if (!is.null(input$selectIndiv) && input$selectIndiv != "ALL") 
      haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv) 
    
    haplo.filter <- haplo.filter %>% rename("Individual ID"=id)
    
    DT::datatable(
      haplo.filter, options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )
    )
  })
  

  output$haploFreqTbl <- DT::renderDataTable({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()   
  
    DT::datatable(
      haplo.freqTbl() %>% mutate(obs.freq=round(obs.freq,3), expected.freq=round(expected.freq,3)), options = list(
          lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
          pageLength = 15
        )
    )
  })

  output$haploSummary <- DT::renderDataTable({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv))
      return ()   
    
    DT::datatable(
      haplo.summaryTbl() %>%
        rename("Individual ID"=id), options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15
      )
    )
  })
  
  
})