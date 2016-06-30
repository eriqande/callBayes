library(shiny)
library("ggplot2")
library("plyr")
library("dplyr")
library("tidyr")
library("DT")
library("grid")
library("scales")

# haplo.cutoff <- haplo.sum %>%
#   group_by(locus, id) %>%
#   summarise(hapl.three.pl.st = ifelse(length(depth) > 2, 0, 0),
#             hapl.three.pl.end = ifelse(length(depth) > 2, sort(depth, decr=T)[3]-1, 0),
#             hapl.one.st = ifelse(sum(depth==max(depth))==1 && length(depth) > 1, sort(depth, decr=T)[2], 0),
#             hapl.one.end = ifelse(sum(depth==max(depth))==1, max(depth),0))


shinyServer(function(input, output, session) {

  dirFiles <- list.files()
  rds.file <- grep(".rds", dirFiles)

  if(length(rds.file)>0){
  select.file.tem <- dirFiles[rds.file[1]]
  updateSelectInput(session, "selectDB", selected=select.file.tem, choices=dirFiles[rds.file])
  haplo.sum<- readRDS(select.file.tem)  %>% mutate(id = as.character(id))
  colnames(haplo.sum) <- c("group", "id", "locus", "haplo", "depth", "logP.call", "logP.miscall", "pos", "allele.balance","rank")
  }
  else {
    haplo.sum <- NULL
  }


  #makeReactiveBinding("haplo.sum")

  update.Haplo.file <- reactive({
    if(input$selectDB == "" || is.null(input$selectDB) || !file.exists(input$selectDB)) return()
    #cat(file=stderr(), "select DB_", input$selectDB, "_----\n")
    readRDS(input$selectDB)  %>% mutate(id = as.character(id))

  })


  ranges <- reactiveValues(y = NULL, x = NULL)
  rangesH <- reactiveValues(y = NULL)
  locusPg <- reactiveValues(l = NULL, width=NULL)
  indivPg <- reactiveValues(i = NULL, width=NULL)
  groupPg <- reactiveValues(g = NULL, width =1)
  hapPg <- reactiveValues(width=NULL)

  filterParam <- reactiveValues(minRead = 1, minAllele = 0.2)
  panelParam <- reactiveValues(n.locus = NULL,
                               n.indiv = NULL,
                               locus.label.tbl = NULL,
                               locus.label = NULL,
                               locus.label.bare =NULL,
                               indiv.label.tbl = NULL,
                               indiv.label = NULL,
                               indiv.label.bare = NULL,
                               is.reject=NULL,
                               n.group = 0,
                               group.label = NULL,
                               group.label.tbl = NULL,
                               group.label.bare = NULL)


  observeEvent(input$selectDB, {
    cat(file=stderr(), "select DB_", input$selectDB, "_","----\n")
    if(is.null(haplo.sum) || input$selectDB == "" || is.null(input$selectDB)) return()

    haplo.sum <- update.Haplo.file()

    if (!"group" %in% colnames(haplo.sum)) haplo.sum <- cbind.data.frame("group"="unlabel",haplo.sum, stringsAsFactors=F) %>% tbl_df
    cat(file=stderr(), "preview_", head(haplo.sum,1) %>% unlist(), "_----\n")
    panelParam$n.locus <- length(unique(haplo.sum$locus))
    panelParam$n.indiv <- length(unique(haplo.sum$id))

    locus.sorted <- sort(unique(haplo.sum$locus))
    panelParam$locus.label.tbl <-  data.frame(locus =locus.sorted, stringsAsFactors = F) %>% tbl_df()
    panelParam$locus.label <- c("ALL",locus.sorted)
    panelParam$locus.label.bare <- locus.sorted

    indiv.sorted <- sort(unique(haplo.sum$id))
    panelParam$indiv.label.tbl <-  data.frame(id =indiv.sorted, stringsAsFactors = F) %>% tbl_df()
    panelParam$indiv.label <- c("ALL",indiv.sorted)
    panelParam$indiv.label.bare <- indiv.sorted

    group.sorted <- sort(unique(haplo.sum$group))
    panelParam$group.label.tbl <-  data.frame(id =group.sorted, stringsAsFactors = F) %>% tbl_df()
    panelParam$group.label <- c("ALL",group.sorted)
    panelParam$group.label.bare <- group.sorted

    panelParam$is.reject <- rep(0, panelParam$n.locus)

    updateSelectInput(session, "selectLocus", selected="ALL", choices=panelParam$locus.label)
    updateSelectInput(session, "selectIndiv", selected="ALL", choices=panelParam$indiv.label)
    updateSelectInput(session, "selectGroup", selected="ALL", choices=panelParam$group.label)


    end.indx<-min(15,panelParam$n.locus)
    locusPg$l <- panelParam$locus.label.bare[1:end.indx]
    rangesH$y <- c(0, length(locusPg$l)+1)

    end.indx<-min(15,panelParam$n.indiv)
    indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
    ranges$y <- c(0, length(indivPg$i)+1)

    Filter.haplo.sum()

  },priority = -3)

  ## updating Locus and individidual choice at the start of the session:
  #updateSelectInput(session, "selectLocus", selected="ALL", choices=locus.label)
  #updateSelectInput(session, "selectIndiv", selected="ALL", choices=indiv.label)

  # reacting to the locus & Indiv's previous and next button
  observeEvent(input$locusBack, {
    indx <- isolate(which(panelParam$locus.label==input$selectLocus))
    label <- ifelse(indx > 1, panelParam$locus.label[indx-1], panelParam$locus.label[indx])
    updateSelectInput(session, "selectLocus", selected=label)
  })
  observeEvent(input$locusFor, {
    indx <- isolate(which(panelParam$locus.label==input$selectLocus))
    label <- ifelse(indx < length(panelParam$locus.label), panelParam$locus.label[indx+1], panelParam$locus.label[indx])
    updateSelectInput(session, "selectLocus", selected=label)
  })
  observeEvent(input$indivBack, {
    indx <- isolate(which(panelParam$indiv.label==input$selectIndiv))
    label <- ifelse(indx > 1, panelParam$indiv.label[indx-1], panelParam$indiv.label[indx])
    updateSelectInput(session, "selectIndiv", selected=label)
  })
  observeEvent(input$indivFor, {
    indx <- isolate(which(panelParam$indiv.label==input$selectIndiv))
    label <- ifelse(indx < length(panelParam$indiv.label), panelParam$indiv.label[indx+1], panelParam$indiv.label[indx])
    updateSelectInput(session, "selectIndiv", selected=label)
  })

  observeEvent(input$groupBack, {
    indx <- isolate(which(panelParam$group.label==input$selectGroup))
    label <- ifelse(indx > 1, panelParam$group.label[indx-1], panelParam$group.label[indx])
    updateSelectInput(session, "selectGroup", selected=label)
  })
  observeEvent(input$groupFor, {
    indx <- isolate(which(panelParam$group.label==input$selectGroup))
    label <- ifelse(indx < length(panelParam$group.label), panelParam$group.label[indx+1], panelParam$group.label[indx])
    updateSelectInput(session, "selectGroup", selected=label)
  })


  # reacting to the filter update button
  observeEvent(input$updateFilter, {
    filterParam$minRead <- input$coverageMin
    filterParam$minAllele <- input$minAlleleRatio
    Filter.haplo.sum()
  })

  observeEvent(input$selectGroup, {
    if(is.null(haplo.sum) || input$selectDB == "" || is.null(input$selectDB)) return()
    indx <- isolate(which(panelParam$group.label.bare==input$selectGroup))
    haplo.sum <- update.Haplo.file()
    if (input$selectGroup != "ALL") haplo.sum <- haplo.sum %>% filter(group == input$selectGroup)

    panelParam$n.indiv <- length(unique(haplo.sum$id))

    indiv.sorted <- sort(unique(haplo.sum$id))
    panelParam$indiv.label.tbl <-  data.frame(id =indiv.sorted, stringsAsFactors = F) %>% tbl_df()
    panelParam$indiv.label <- c("ALL",indiv.sorted)
    panelParam$indiv.label.bare <- indiv.sorted
    updateSelectInput(session, "selectIndiv", selected="ALL", choices=panelParam$indiv.label)

    end.indx<-min(15,panelParam$n.indiv)
    indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
    ranges$y <- c(0, length(indivPg$i)+1)

    Filter.haplo.sum()
  }, priority = -2)


  observeEvent(input$selectLocus,{
    indx <- isolate(which(panelParam$locus.label.bare==input$selectLocus))

    output$locusSelect <- renderText({input$selectLocus})
    output$locusSelect1 <- renderText({input$selectLocus})
    if(input$selectLocus != "ALL"){
      output$maxlocusPage <- renderText({"1"})
      updateNumericInput(session, "locusPage", value=1, max=1)
      locusPg$l <- input$selectLocus
      rangesH$y <- c(0, 2)
      output$locusAcceptStatus <- renderText({ifelse(panelParam$is.reject[indx]==0,"Accept","Reject")})
    }
    else {

      output$locusAcceptStatus <- renderText({"NA"})
      output$maxlocusPage <- renderText({paste0(ceiling(as.numeric(panelParam$n.locus)/
                                                          as.numeric(input$locusPerDisplay)))})
      updateNumericInput(session, "locusPage", value=1, max=ceiling(panelParam$n.locus/15))
      updateSelectInput(session, "locusPerDisplay", selected = 15)
      end.indx<-min(15,panelParam$n.locus)
      locusPg$l <- panelParam$locus.label.bare[1:end.indx]
      rangesH$y <- c(0, end.indx+1)
    }
    Filter.haplo.sum()
  })

  observeEvent(input$locusPerDisplay,{
    if (is.null(panelParam$n.locus)) return()

    if(input$selectLocus != "ALL"){
      output$maxlocusPage <- renderText({"1"})
      updateNumericInput(session, "locusPage", value=1, max=1)
      locusPg$l <- input$selectLocus
      rangesH$y <- c(0, 2)
    }
    else {
      if (input$locusPerDisplay == 100) {
        output$maxlocusPage <- renderText({"1"})
        updateNumericInput(session, "locusPage", value=1, max=1)
        locusPg$l <- panelParam$locus.label.bare
        rangesH$y <- c(0, length(locusPg$l)+1)
      }
      else{
        output$maxlocusPage <- renderText({paste0(ceiling(as.numeric(panelParam$n.locus)/
                                                            as.numeric(input$locusPerDisplay)))})
        #cat(file=stderr(), "haha_", as.numeric(panelParam$n.locus)/as.numeric(input$locusPerDisplay), "_----\n")
        updateNumericInput(session, "locusPage", max=ceiling(as.numeric(panelParam$n.locus)
                                                             /as.numeric(input$locusPerDisplay)))

        end.indx<-min(as.numeric(panelParam$n.locus) ,as.numeric(input$locusPerDisplay))
        locusPg$l <- panelParam$locus.label.bare[1:end.indx]
        rangesH$y <- c(0, length(locusPg$l)+1)
      }
    }

  })

  observeEvent(input$updateLocusSizeDisplay, {
    if(input$selectLocus == "ALL"){
      if (input$locusPerDisplay == 100) {
        locusPg$l <- panelParam$locus.label.bare
        rangesH$y <- c(0, length(locusPg$l)+1)

      }
      else {
        pg <- min(as.numeric(input$locusPage), ceiling(as.numeric(panelParam$n.locus)/
                                                         as.numeric(input$locusPerDisplay)))
        start.indx <- (as.numeric(input$locusPerDisplay)*(pg-1))+1
        end.indx<-min(as.numeric(panelParam$n.locus) ,as.numeric(input$locusPerDisplay)*pg)
        locusPg$l <- panelParam$locus.label.bare[start.indx:end.indx]
        rangesH$y <- c(0, length(locusPg$l)+1)

      }
    }
  })


  observeEvent(input$selectIndiv,{
    output$indivSelect <- renderText({input$selectIndiv})
    if(input$selectIndiv != "ALL"){
      output$maxIndivPage <- renderText({"1"})
      updateNumericInput(session, "indivPage", value=1, max=1)
      indivPg$i <- input$selectIndiv
      ranges$y <- c(0, 2)
    }
    else {
      output$maxIndivPage <- renderText({paste0(ceiling(as.numeric(panelParam$n.indiv)/
                                                          as.numeric(input$indivPerDisplay)))})
      updateNumericInput(session, "indivPage", value=1, max=ceiling(panelParam$n.indiv/15))
      updateSelectInput(session, "indivPerDisplay", selected = 15)
      end.indx<-min(15,panelParam$n.indiv)
      indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
      ranges$y <- c(0, end.indx+1)
    }
    Filter.haplo.sum()
  })

  observeEvent(input$acceptLocus,{
    if(input$selectLocus != "ALL"){
      indx <- isolate(which(panelParam$locus.label.bare==input$selectLocus))
      panelParam$is.reject[indx] <- 0
      output$locusAcceptStatus <- renderText({ifelse(panelParam$is.reject[indx]==0,"Accept","Reject")})
    }
  })

  observeEvent(input$rejectLocus,{
    if(input$selectLocus != "ALL"){
      indx <- isolate(which(panelParam$locus.label.bare==input$selectLocus))
      panelParam$is.reject[indx] <- 1
      output$locusAcceptStatus <- renderText({ifelse(panelParam$is.reject[indx]==0,"Accept","Reject")})
    }
  })

  observeEvent(input$indivPerDisplay,{
    if (is.null(panelParam$n.indiv)) return()

    if(input$selectIndiv != "ALL"){
      output$maxIndivPage <- renderText({"1"})
      updateNumericInput(session, "indivPage", value=1, max=1)
      indivPg$i <- input$selectIndiv
      ranges$y <- c(0, 2)
    }
    else {
      if (input$indivPerDisplay == 100) {
        output$maxIndivPage <- renderText({"1"})
        updateNumericInput(session, "indivPage", value=1, max=1)
        indivPg$i <- panelParam$indiv.label.bare
        ranges$y <- c(0, length(indivPg$i)+1)
      }
      else{
        output$maxIndivPage <- renderText({paste0(ceiling(as.numeric(panelParam$n.indiv)/
                                                            as.numeric(input$indivPerDisplay)))})
        updateNumericInput(session, "indivPage", max=ceiling(as.numeric(panelParam$n.indiv)
                                                             /as.numeric(input$indivPerDisplay)))

        end.indx<-min(as.numeric(panelParam$n.indiv) ,as.numeric(input$indivPerDisplay))
        indivPg$i <- panelParam$indiv.label.bare[1:end.indx]
        ranges$y <- c(0, length(indivPg$i)+1)
      }
    }

  })

  observeEvent(input$updateIndivSizeDisplay, {
    if(input$selectIndiv == "ALL"){
      if (input$indivPerDisplay == 100) {
        indivPg$i <- panelParam$indiv.label.bare
        ranges$y <- c(0, length(indivPg$i)+1)

      }
      else {
        pg <- min(as.numeric(input$indivPage), ceiling(as.numeric(panelParam$n.indiv)/
                                                         as.numeric(input$indivPerDisplay)))
        start.indx <- (as.numeric(input$indivPerDisplay)*(pg-1))+1
        end.indx<-min(as.numeric(panelParam$n.indiv) ,as.numeric(input$indivPerDisplay)*pg)
        indivPg$i <- panelParam$indiv.label.bare[start.indx:end.indx]
        ranges$y <- c(0, length(indivPg$i)+1)

      }
    }
  })


  haplo.summaryTbl <- reactive({
    haplo.sum <- update.Haplo.file()
    if(is.null(haplo.sum)) return ()
    haplo.filter <- haplo.sum %>%
      filter(depth > filterParam$minRead, rank <= 2, allele.balance >= filterParam$minAllele)

    if (input$selectGroup != "ALL") haplo.filter <- haplo.filter %>% filter(group == input$selectGroup)
    if (input$selectLocus != "ALL") haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus)
    if (input$selectIndiv != "ALL") haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)

    haplo.filter <- haplo.filter %>%
      group_by(locus, id, group) %>%
      arrange(-depth) %>%
      summarise(haplotype.1 = ifelse(length(depth)==1, haplo[1], sort(haplo)[1]),
                haplotype.2 = ifelse(length(depth)==1, haplo[1], sort(haplo)[2]),
                read.depth.1 = depth[1],
                read.depth.2 = ifelse(length(depth)==1, depth[1], depth[2]))
  })

  haplo.freqTbl <- reactive({

    if(is.null(haplo.summaryTbl())) {return()}

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


  Filter.haplo.sum <- reactive({
    haplo.sum <- update.Haplo.file()

    if(is.null(haplo.sum)) return ()

    haplo.filter <- haplo.sum %>%
      filter(depth > filterParam$minRead, allele.balance >= filterParam$minAllele)

    if(input$topTwo)
      haplo.filter <- haplo.filter %>% filter(rank <= 2)

    if (input$selectGroup != "ALL") haplo.filter <- haplo.filter %>% filter(group == input$selectGroup)
    if (input$selectLocus != "ALL") haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus)
    if (input$selectIndiv != "ALL") haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)

    groupPg$width <- dim(haplo.filter)[1]
    haplo.filter
  })

  Get.tbl.by.locus <- reactive({
    if (is.null(Filter.haplo.sum())) return()
    haplo.ct <- Filter.haplo.sum() %>%
      group_by(locus, id, group) %>%
      summarise(tot.hapl = n(), tot.depth = sum(depth))

  })

  Get.tbl.by.id <- reactive({
    if (is.null(Filter.haplo.sum())) return()
    haplo.ct <- Filter.haplo.sum() %>%
      group_by(id, locus) %>%
      summarise(tot.depth = sum(depth))
  })



  # BY LOCUS PANEL::

  output$haplDensityPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()

    if(is.null(Get.tbl.by.locus())) return()
    if (dim(panelParam$locus.label.tbl)[1]==0) return()

    haplo.tot.tbl <- Get.tbl.by.locus() %>%
      group_by(locus, tot.hapl) %>%
      summarise(ct = n()) %>%
      ungroup() %>%
      group_by(locus) %>%
      mutate(frac = ct/sum(ct))

    #cat(file=stderr(), "is it updating_", unlist(haplo.tot.tbl[1,]), "_----\n")



    uniqH.perI.tbl <- right_join(haplo.tot.tbl, panelParam$locus.label.tbl, by="locus")
    if(is.null(uniqH.perI.tbl)) return()
    uniqH.perI.tbl[is.na(uniqH.perI.tbl)]<- 0

    #cat(file=stderr(), "is it moving :_ _----\n")


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
      ylim(locusPg$l)+
      coord_cartesian(ylim=rangesH$y)
  },height = function() ifelse(groupPg$width==0,0,
                               max(ifelse(input$selectLocus=="ALL",
                                          ifelse(input$locusPerDisplay==100,
                                                 9*length(panelParam$locus.label),
                                                 9*as.numeric(input$locusPerDisplay)),
                                          1),250)))

  output$numHapPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || is.null(locusPg$l))
      return ()
    if(is.null(haplo.summaryTbl())) {return()}
    if(dim(panelParam$locus.label.tbl)[1]==0) return()


    frac.calleable <- haplo.summaryTbl() %>% group_by(locus) %>% summarise(n=length(unique(c(haplotype.1,haplotype.2))))

    frac.calleable <- right_join(frac.calleable, panelParam$locus.label.tbl, by="locus")
    if(is.null(frac.calleable)) return()

    frac.calleable[is.na(frac.calleable)]<- 0

    if (input$selectLocus != "ALL") {
      frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    }


    ggplot(frac.calleable, aes(x=n, y=locus))+
      geom_point()+
      xlab("num of unique haplotypes")+
      ylab("")+
      scale_size_continuous(guide=FALSE)+#"fraction")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 2, 0, 0), "mm"))+
      ylim(locusPg$l)+
      coord_cartesian(ylim=rangesH$y)
    #scale_x_discrete(limits=c(-1, max(frac.calleable$n)+1)) #breaks= pretty_breaks()
  },height = function() ifelse(groupPg$width==0,0,max(ifelse(input$selectLocus=="ALL",
                                                             ifelse(input$locusPerDisplay==100,
                                                                    9*length(panelParam$locus.label),
                                                                    9*as.numeric(input$locusPerDisplay)),
                                                             1),250)))

  output$fracIndivPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || is.null(locusPg$l))
      return ()
    if(is.null(haplo.summaryTbl())) {return()}
    if(dim(panelParam$locus.label.tbl)[1]==0) return()

    nIndiv <- ifelse(input$selectIndiv == "ALL", panelParam$n.indiv, 1)

    frac.calleable <- haplo.summaryTbl() %>% group_by(locus) %>% summarise(f=n()/nIndiv)
    frac.calleable <- right_join(frac.calleable, panelParam$locus.label.tbl, by="locus")
    if(is.null(frac.calleable)) return()

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
      ylim(locusPg$l)+
      coord_cartesian(ylim=rangesH$y)+
      xlim(c(0,1))
  },height = function() ifelse(groupPg$width==0,0,
                               max(ifelse(input$selectLocus=="ALL",
                                          ifelse(input$locusPerDisplay==100,
                                                 9*length(panelParam$locus.label),
                                                 9*as.numeric(input$locusPerDisplay)),
                                          1),250)))

  output$readDepthPerLocus <- renderPlot({

    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()

    if(is.null(Get.tbl.by.locus())) return()
    if (dim(panelParam$locus.label.tbl)[1]==0) return()

    readDepth.perI.tbl <- right_join(Get.tbl.by.locus(), panelParam$locus.label.tbl, by="locus")
    if(is.null(readDepth.perI.tbl)) return()
    readDepth.perI.tbl[is.na(readDepth.perI.tbl)]<- 0

    if (input$selectLocus != "ALL") {
      readDepth.perI.tbl <- readDepth.perI.tbl %>% filter(locus == input$selectLocus)
    }

    readDepth.perI.tbl <- readDepth.perI.tbl %>% group_by(locus) %>% mutate(mean.depth = mean(tot.depth))

    ggplot(readDepth.perI.tbl, aes(x=locus, y=tot.depth)) +
      xlab("")+
      ylab("read depth per individual")+
      geom_violin()+
      geom_point(aes(x=locus, y=mean.depth), cex=3, pch=3)+
      #geom_point()+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      scale_y_log10()+
      xlim(locusPg$l)+
      coord_flip(xlim=rangesH$y)
  },height = function() ifelse(groupPg$width==0,0,max(ifelse(input$selectLocus=="ALL",
                                                             ifelse(input$locusPerDisplay==100,
                                                                    9*length(panelParam$locus.label),
                                                                    9*as.numeric(input$locusPerDisplay)),
                                                             1),250)))



  ## BY INDIVIDUAL PANEL::

  output$AlleleRatioByIndiv <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()

    haplo.sum <- update.Haplo.file()

    haplo.filter <- haplo.sum %>%
      filter(depth > filterParam$minRead, allele.balance >= filterParam$minAllele)

    if(input$topTwo)
      haplo.filter <- haplo.filter %>% filter(rank <= 2)

    if (input$selectLocus != "ALL")
      haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus)

    if (dim(haplo.filter)[1]==0 ) return ()
    if (dim(panelParam$indiv.label.tbl)[1]==0) return()

    haplo.filter <- haplo.filter %>%
      group_by(locus, id, group) %>%
      summarise(depth.ratio = ifelse(length(depth)==1, 0, min(allele.balance)),
                depth.first = max(depth))

    haplo.filter <- right_join(haplo.filter, panelParam$indiv.label.tbl, by="id")
    if(is.null(haplo.filter)) return()
    haplo.filter[is.na(haplo.filter)]<- 0

    #if (input$selectIndiv != "ALL") {
    #  haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #}

    ggplot(data=haplo.filter, aes(x=depth.ratio, y=id,size=log(depth.first, 10), color=group ))+
      geom_point(alpha=0.4)+
      scale_size_continuous(guide=FALSE)+#"Read Depth of the most common haplotype (log 10)")+
      scale_color_discrete(guide=FALSE, drop=T, limits=levels(panelParam$group.label.bare))+
      theme_bw()+
      ylab("individual ID")+
      xlab ("depth ratio of the second common haplotype : first common haplotype")+
      theme(legend.position="bottom")+
      xlim(c(0,1))+
      ylim(indivPg$i)+
      coord_cartesian(ylim=ranges$y)+
      geom_vline(xintercept=filterParam$minAllele, linetype="dashed", color = "red")
  },height = function()  ifelse(groupPg$width==0,0,
                                max(ifelse(input$selectIndiv=="ALL",
                                           ifelse(input$indivPerDisplay==100,
                                                  9*length(panelParam$indiv.label),
                                                  9*as.numeric(input$indivPerDisplay)),
                                           1),250)))



  output$numUniqHapByIndiv <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || is.null(indivPg$i))
      return ()
    if(dim(panelParam$locus.label.tbl)[1]==0) return()

    filter.haplo <- Filter.haplo.sum() %>%
      group_by(id, locus) %>%
      summarise(n.hap.locus=n()) %>%
      ungroup() %>%
      group_by(id, n.hap.locus) %>%
      summarise(n.locus=n())


    tot.hap.per.indiv <- right_join(filter.haplo, panelParam$indiv.label.tbl, by="id")
    if(is.null(tot.hap.per.indiv)) return()

    tot.hap.per.indiv[is.na(tot.hap.per.indiv)]<- 0


    ggplot( tot.hap.per.indiv, aes(x=n.hap.locus, y=id, size=n.locus))+
      geom_point()+
      xlab("num of haplotypes per locus")+
      ylab("")+
      scale_size_continuous(guide=FALSE)+#"fraction")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.margin = unit(0, 'mm'),
        plot.margin = unit(c(0, 2, 0, 0), "mm"))+
      ylim(indivPg$i)+
      coord_cartesian(ylim=ranges$y)
    #scale_x_discrete(limits=c(-1, max(frac.calleable$n)+1)) #breaks= pretty_breaks()
  },height = function()  ifelse(groupPg$width==0,0,
    max(ifelse(input$selectIndiv=="ALL",
      ifelse(input$indivPerDisplay==100,
        9*length(panelParam$indiv.label),
        9*as.numeric(input$indivPerDisplay)),
      1),250)))



  output$fracHaploPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()
    if(is.null(haplo.summaryTbl())) {return()}

    nLocus <- ifelse(input$selectLocus == "ALL", panelParam$n.locus, 1)

    haplo.filter <- haplo.summaryTbl() %>% group_by(id) %>% summarise(f=n()/nLocus)
    if(dim(panelParam$indiv.label.tbl)[1]==0) return()
    haplo.filter <- right_join(haplo.filter, panelParam$indiv.label.tbl, by="id")
    if(is.null(haplo.filter)) return()

    haplo.filter[is.na(haplo.filter)]<- 0

    #if (input$selectIndiv != "ALL") {
    #  haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #}

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
      ylim(indivPg$i)+
      coord_cartesian(ylim=ranges$y)+
      xlim(c(0,1))
  },height = function()  ifelse(groupPg$width==0,0,
                                max(ifelse(input$selectIndiv=="ALL",
                                           ifelse(input$indivPerDisplay==100,
                                                  9*length(panelParam$indiv.label),
                                                  9*as.numeric(input$indivPerDisplay)),
                                           1),250)))

  # output$meanReadDepthByIndiv <- renderPlot({
  #   if (is.null(input$selectLocus) || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
  #     return ()
  #
  #   if(is.null(Get.tbl.by.id())) return()
  #
  #   haplo.filter <- Get.tbl.by.id() %>%
  #     ungroup()%>%
  #     group_by(id) %>%
  #     summarise(mean.depth = mean(tot.depth))
  #
  #   if(dim(panelParam$indiv.label.tbl)[1]==0) return()
  #   haplo.filter <- right_join(haplo.filter, panelParam$indiv.label.tbl, by="id")
  #   if(is.null(haplo.filter)) return()
  #   haplo.filter[is.na(haplo.filter)]<- 0.0#0.0001 if turn on log 10 scale
  #
  #   #if (input$selectIndiv != "ALL") {
  #   #  haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
  #   #}
  #
  #   ggplot(haplo.filter, aes(y=id, x=mean.depth)) +
  #     geom_point()+
  #     ylab("")+
  #     xlab("mean locus read depth ")+
  #     theme_bw()+
  #     theme(axis.text.y=element_blank(),
  #           axis.ticks.y=element_blank(),
  #           panel.margin = unit(0, 'mm'))+
  #     #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  #     #scale_x_log10()+
  #     ylim(indivPg$i)+
  #     coord_cartesian(ylim=ranges$y)
  #
  # },height = function() ifelse(groupPg$width==0,0,
  #                              max(ifelse(input$selectIndiv=="ALL",
  #                                         ifelse(input$indivPerDisplay==100,
  #                                                9*length(panelParam$indiv.label),
  #                                                9*as.numeric(input$indivPerDisplay)),
  #                                         1),250)))

  output$readDepthByIndiv <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()
    if (is.null(Filter.haplo.sum())) return()
    if(dim(panelParam$indiv.label.tbl)[1]==0) return()



    haplo.filter <- right_join( Filter.haplo.sum(), panelParam$indiv.label.tbl, by="id")
    if(is.null(haplo.filter)) return()
    haplo.filter[is.na(haplo.filter)]<- 0

    haplo.filter <- haplo.filter %>% group_by(id) %>% mutate(mean.depth=mean(depth))

    #  if (input$selectIndiv != "ALL") {
    #    haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #  }

    ggplot(haplo.filter, aes(x=id, y=depth, group=id)) +
      xlab("")+
      ylab("haplotype read depth")+
      geom_violin()+
      geom_point(aes(x=id, y=mean.depth),pch=3, cex=3)+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'))+
      #plot.margin = unit(c(0, 0, 0, 0), "mm"))+
      #scale_y_log10()+
      xlim(indivPg$i)+
      coord_flip(xlim=ranges$y)

  },height = function() ifelse(groupPg$width==0,0,
                               max(ifelse(input$selectIndiv=="ALL",
                                          ifelse(input$indivPerDisplay==100,
                                                 9*length(panelParam$indiv.label),
                                                 9*as.numeric(input$indivPerDisplay)),
                                          1),250)))

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
      ranges$y <- c(0, length(indivPg$i)+1)
    }
  })

  #   observeEvent(input$plot1_dblclick, {
  #     brush <- input$plot1_brush
  #     if (!is.null(brush)) {
  #       ranges$y <- c(brush$ymin, brush$ymax)
  #       ranges$x <- c(brush$xmin, brush$xmax)
  #     } else {
  #       ranges$y <- NULL
  #       ranges$x <- NULL
  #     }
  #   })

  observeEvent(input$plotH_dblclick, {
    brush <- input$plotH_brush
    if (!is.null(brush)) {
      rangesH$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesH$y <- c(0, length(locusPg$l)+1)
    }
  })


  # by-group distribution panel

  output$nIndivByGroupPlot <- renderPlot({

    if (is.null(input$selectLocus) || is.null(input$selectIndiv) ||  input$selectDB == "" || is.null(locusPg$l) )
      return ()

    filter.tbl <- Filter.haplo.sum() %>% group_by(group) %>% summarise(n.indiv=length(unique(id)))

    ggplot(filter.tbl, aes(x=n.indiv, y=group, color=group))+
      geom_point()+
      xlab("num of indiv")+
      ylab("")+
      scale_color_discrete(guide=FALSE)+
      theme_bw()+
      theme(#axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.margin = unit(0, 'mm'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))
  },height = function(){ifelse(groupPg$width==0,0,
                               ifelse(input$selectGroup=="ALL",300,100)) })

  output$fIndivByGroupPlot <- renderPlot({

    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || input$selectDB == "" )
      return ()


    all.indiv <- update.Haplo.file() %>% group_by(group, locus) %>% summarise(nIndiv=ifelse( input$selectIndiv!="ALL",1,length(unique(id))) )
    filter.indiv <- Filter.haplo.sum() %>% group_by(group, locus) %>% summarise(fIndiv=length(unique(id)))
    frac.calleable <- left_join(filter.indiv, all.indiv, by=c("group", "locus")) %>% mutate(f=fIndiv/nIndiv)

    mean.f.tbl <- frac.calleable %>% group_by(group) %>% summarise(mean.f = mean(f, na.rm =T))
    #     if(is.null(frac.calleable)) return()
    #     frac.calleable[is.na(frac.calleable)]<- 0
    #
    #     if (input$selectLocus != "ALL") {
    #       frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    #     }

    ggplot(frac.calleable, aes(x=f, y=group, color=group))+
      geom_point(alpha=0.5)+
      geom_point(data=mean.f.tbl, aes(y=group, x = mean.f), color = "black", pch=3, cex=3)+
      xlab("fraction of indiv w/ calleable haplotype")+
      ylab("")+
      scale_color_discrete(guide=FALSE)+#"fraction")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 0, 0, 0), "mm"))
  },height = function(){ifelse(groupPg$width==0,0,
                               ifelse(input$selectGroup=="ALL",300,100)) })

  output$nLociByGroupPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || input$selectDB == "" || is.null(locusPg$l) )
      return ()

    filter.tbl <- Filter.haplo.sum() %>% group_by(group) %>% summarise(n.locus=length(unique(locus)))

    ggplot(filter.tbl, aes(x=n.locus, y=group, color=group))+
      geom_point()+
      xlab("num of locus")+
      ylab("")+
      scale_color_discrete(guide=FALSE)+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 0, 0, 0), "mm"))
  },height = function(){ifelse(groupPg$width==0,0,
                               ifelse(input$selectGroup=="ALL",300,100)) })

  output$fLociByGroupPlot <- renderPlot({
    if (is.null(input$selectLocus) || is.null(input$selectIndiv) || input$selectDB == "" || is.null(locusPg$l) )
      return ()

    all.locus <- update.Haplo.file() %>% group_by(group, id) %>% summarise(nLocus=ifelse( input$selectLocus!="ALL",1,length(unique(locus))) )
    filter.locus <- Filter.haplo.sum() %>% group_by(group, id) %>% summarise(fLocus=length(unique(locus)))
    frac.calleable <- left_join(filter.locus, all.locus, by=c("group", "id")) %>% mutate(f=fLocus/nLocus)

    mean.f.tbl <- frac.calleable %>% group_by(group) %>% summarise(mean.f = mean(f, na.rm =T))
    #     if(is.null(frac.calleable)) return()
    #     frac.calleable[is.na(frac.calleable)]<- 0
    #
    #     if (input$selectLocus != "ALL") {
    #       frac.calleable <- frac.calleable %>% filter(locus == input$selectLocus)
    #     }

    ggplot(frac.calleable, aes(x=f, y=group, color=group))+
      geom_point(alpha=0.5)+
      geom_point(data=mean.f.tbl, aes(y=group, x = mean.f), color = "black", pch=3, cex=3)+
      xlab("fraction of loci w/ calleable haplotype")+
      ylab("")+
      scale_color_discrete(guide=FALSE)+#"fraction")+
      theme_bw()+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.margin = unit(0, 'mm'),
            plot.margin = unit(c(0, 0, 0, 0), "mm"))
  },height = function(){ifelse(groupPg$width==0,0,
                               ifelse(input$selectGroup=="ALL",300,100)) })

  #ABOUT HAPLOTYPE distribution panel
  output$hapSeq <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()

    haplo.filter <- Filter.haplo.sum()
    position <- strsplit(haplo.filter$pos[1],",") %>% unlist

    haplo.profile.frac <- haplo.filter %>%
      group_by(haplo) %>%
      summarise(n=n()) %>%
      ungroup() %>%
      mutate(frac = n/sum(n))

    if(nrow(haplo.profile.frac)==0) return()

    haplo.split.profile <- sapply(1:nrow(haplo.profile.frac), function(i) {
      char.split <- strsplit(haplo.profile.frac[i,]$haplo, "")
      #cat(file=stderr(), "character split_", unlist(char.split), "_----\n")
      n.char <- length(char.split[[1]])
      sapply(1:n.char, function(j) c(i, j, char.split[[1]][j], haplo.profile.frac[i,]$frac))
    }) %>%
      matrix(., ncol=4, byrow=T) %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      tbl_df()

    colnames(haplo.split.profile) <- c("group", "pos", "seq", "frac")
    haplo.split.profile <- haplo.split.profile %>% mutate(pos=as.numeric(position[as.numeric(pos)]),
                                                          frac=as.numeric(frac),
                                                          group=as.numeric(group))


    g <- ggplot(haplo.split.profile, aes(x=pos, y=seq, group=group, size=frac, color=factor(group))) +
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
  },height = function(){ifelse(groupPg$width==0,0, ifelse(input$selectLocus=="ALL",0,400)) })

  output$histHap <- renderPlot({

    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()

    #     haplo.sum <- update.Haplo.file()
    #     if(is.null(haplo.sum)) return ()
    #
    #     nIndiv <- ifelse(input$selectIndiv == "ALL", panelParam$n.indiv, 1)
    #
    #     haplo.filter <- haplo.sum %>%
    #       filter(depth > filterParam$minRead, locus == input$selectLocus, allele.balance >= filterParam$minAllele)
    #     if(input$topTwo)
    #       haplo.filter <- haplo.filter %>% filter(rank <= 2)
    #     if (input$selectIndiv != "ALL")
    #       haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    #     haplo.filter <- haplo.filter %>% group_by(haplo) %>% summarise(f=n()/nIndiv)
    if(is.null(haplo.summaryTbl())) {return()}

    obs.freq.tbl<-  haplo.summaryTbl() %>%
      ungroup() %>%
      group_by(locus) %>%
      mutate(tot.haplo = n()) %>%
      group_by(locus, haplotype.1, haplotype.2) %>%
      summarise(obs.freq=n()/tot.haplo[1])

    allelic.freq.tbl <- gather(obs.freq.tbl, whichHap,  hap1, 2:3) %>%
      group_by(locus, hap1) %>%
      summarise(f=sum(obs.freq/2))

    if(nrow(allelic.freq.tbl)==0) return()

    ggplot(allelic.freq.tbl, aes(y=hap1, x = f, color=factor(hap1))) +
      geom_point(size=4)+
      xlab("observed freq")+
      ylab("haplotype")+
      theme_bw()+
      scale_color_discrete(guide=FALSE)

  },height = function(){ifelse(groupPg$width==0,0,
                               ifelse(input$selectLocus=="ALL",0,400)) })



  output$PairWiseHap <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()

    haplo.sum <- update.Haplo.file()
    if(is.null(haplo.sum)) return ()


    haplo.filter <- haplo.sum %>%
      filter(depth > filterParam$minRead, locus == input$selectLocus, allele.balance >= filterParam$minAllele)

    #if(input$topTwo)
    haplo.filter <- haplo.filter %>% filter(rank <= 2)
    if (input$selectIndiv != "ALL") haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
    if (input$selectGroup != "ALL") haplo.filter <- haplo.filter %>% filter(group == input$selectGroup)

    if(nrow(haplo.filter)==0) return()


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
      mutate(re.hap1 = sort(c(hap1, hap2))[1],
             re.hap2 = sort(c(hap1, hap2))[2])

    ggplot(haplo.filter, aes(x=hap1, y=hap2, size=n, color=hap1==hap2))+
      geom_point()+
      xlab("haplotype 1")+
      ylab("haplotype 2")+
      geom_point(data=freq.hap, aes(x=re.hap1,y=re.hap2, size=freq*n.hap/2), shape=21, fill=NA, color="black")+
      scale_color_discrete(guide=FALSE)+
      scale_size_continuous(range = c(3,20),guide=FALSE)+
      theme_bw()
  },height = function(){ifelse(groupPg$width==0,0,
                               ifelse(input$selectLocus=="ALL",0,400)) })

  output$hapByGroupPlot <- renderPlot({
    if (is.null(input$selectLocus) || input$selectLocus == "ALL" || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB) || is.null(locusPg$l) )
      return ()

    if(is.null(haplo.summaryTbl())) {return()}

    obs.freq.tbl<-  haplo.summaryTbl() %>%
      ungroup() %>%
      group_by(group, locus) %>%
      mutate(tot.haplo = n()) %>%
      group_by(group, locus, haplotype.1, haplotype.2) %>%
      summarise(obs.freq=n()/tot.haplo[1])

    allelic.freq.tbl <- gather(obs.freq.tbl, whichHap,  hap1, 3:4) %>%
      group_by(group, locus, hap1) %>%
      summarise(f=sum(obs.freq/2))

    ggplot(allelic.freq.tbl, aes(x=group, y = hap1, color=hap1, size=f)) +
      geom_point()+
      xlab("")+
      ylab("")+
      theme_bw()+
      scale_color_discrete(guide=FALSE)+
      scale_size_continuous(guide=FALSE)

  },height = function(){ifelse(groupPg$width==0,0,
                               ifelse(input$selectLocus=="ALL",0,300)) })





  output$downloadData <- downloadHandler(
    filename = function(){
      if(input$selectTbl=="reported indiv haplotype") return ("filtered_haplotype.csv")
      if(input$selectTbl=="observed variants") return ("observed_haplotype.csv")
      if(input$selectTbl=="SNP report") return ("snp_report.csv")
    }
    ,
    content = function(file) {
      if(input$selectTbl=="reported indiv haplotype") {
        if(is.null(haplo.summaryTbl())) return()
        haplo.freq <- haplo.freqTbl() %>% mutate(obs.freq=round(obs.freq,3), expected.freq=round(expected.freq,3))
        haplo.all.tbl <- haplo.summaryTbl() %>% rename("indiv.ID"=id)
        haplo.all <- left_join(haplo.all.tbl, haplo.freq, by=c("locus", "haplotype.1","haplotype.2"))
        haplo.isAccept <- data.frame(locus=panelParam$locus.label.bare, is.reject=panelParam$is.reject,
                                     stringsAsFactors=FALSE)
        haplo.all <- left_join(haplo.all, haplo.isAccept, by=c("locus"))
        write.csv(haplo.all,file)
      }
      if (input$selectTbl ==  "observed variants") {
        haplo.all <- Filter.haplo.sum() %>% rename("indiv.ID"=id)
        write.csv(haplo.all,file)
      }
      if (input$selectTbl=="SNP report") {
        haplo.summaryTable <- haplo.summaryTbl()
        n.base <- nchar(haplo.summaryTable$haplotype.1)
        haplo.all <- haplo.summaryTable[rep(seq(1, nrow(haplo.summaryTable)),n.base),] %>%
          group_by(locus, id, group) %>%
          mutate(snp.id=row_number(),
                 snp=paste0(substr(haplotype.1,snp.id,snp.id),
                            "/",
                            substr(haplotype.2,snp.id,snp.id))) %>%
          select(-haplotype.1, -haplotype.2, -read.depth.1, -read.depth.2)
        write.csv(haplo.all,file)

      }


    }
  )

  observeEvent(input$updateTable,{
    if (is.null(haplo.freqTbl()) || is.null(haplo.summaryTbl())) return()

    if (input$selectTbl=="reported indiv haplotype") {
      haplo.freq <- haplo.freqTbl() %>% mutate(obs.freq=round(obs.freq,3), expected.freq=round(expected.freq,3))
      haplo.all.tbl <- haplo.summaryTbl() %>% rename("Individual ID"=id)
      haplo.all <- left_join(haplo.all.tbl, haplo.freq, by=c("locus", "haplotype.1","haplotype.2"))
      haplo.isAccept <- data.frame(locus=panelParam$locus.label.bare, is.reject=panelParam$is.reject,
                                   stringsAsFactors=FALSE)
      haplo.all <- left_join(haplo.all, haplo.isAccept, by=c("locus"))
    }

    if (input$selectTbl == "observed variants") {
      haplo.all <- Filter.haplo.sum() %>% rename("Individual ID"=id) %>% select(-logP.call,-logP.miscall)
    }

    if (input$selectTbl=="SNP report") {
      haplo.summaryTable <- haplo.summaryTbl()
      n.base <- nchar(haplo.summaryTable$haplotype.1)
      haplo.all <- haplo.summaryTable[rep(seq(1, nrow(haplo.summaryTable)),n.base),] %>%
        group_by(locus, id, group) %>%
        mutate(snp.id=row_number(),
               snp=paste0(substr(haplotype.1,snp.id,snp.id),
                          "/",
                          substr(haplotype.2,snp.id,snp.id))) %>%
        select(-haplotype.1, -haplotype.2, -read.depth.1, -read.depth.2)
    }


    output$haploTbl <- DT::renderDataTable({
      DT::datatable(
        haplo.all, options = list(
          lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
          pageLength = 15
        )
      )
    })

  })

  #   output$haploTbl <- DT::renderDataTable({
  #     if( input$selectDB == "" || is.null(input$selectDB)) return()
  #     haplo.sum <- update.Haplo.file()
  #     if(is.null(haplo.sum)) return ()
  #
  #
  #     haplo.filter <- haplo.sum %>%
  #       filter(depth > filterParam$minRead) %>%
  #       select(id, locus, haplo, depth)
  #
  #     if (!is.null(input$selectLocus) && input$selectLocus != "ALL")
  #       haplo.filter <- haplo.filter %>% filter(locus == input$selectLocus)
  #
  #     if (!is.null(input$selectIndiv) && input$selectIndiv != "ALL")
  #       haplo.filter <- haplo.filter %>% filter(id == input$selectIndiv)
  #
  #     haplo.filter <- haplo.filter %>% rename("Individual ID"=id)
  #
  #     DT::datatable(
  #       haplo.filter, options = list(
  #         lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
  #         pageLength = 15
  #       )
  #     )
  #   })
  #
  #
  #   output$haploFreqTbl <- DT::renderDataTable({
  #     if (is.null(input$selectLocus) || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB))
  #       return ()
  #
  #     if (is.null(haplo.freqTbl()))
  #       return()
  #
  #     DT::datatable(
  #       haplo.freqTbl() %>% mutate(obs.freq=round(obs.freq,3), expected.freq=round(expected.freq,3)), options = list(
  #         lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
  #         pageLength = 15
  #       )
  #     )
  #   })
  #
  #   output$haploSummary <- DT::renderDataTable({
  #     if (is.null(input$selectLocus) || is.null(input$selectIndiv)|| input$selectDB == "" || is.null(input$selectDB))
  #       return ()
  #
  #     if (is.null(haplo.summaryTbl()))
  #       return ()
  #
  #     haplo.freqTbl() %>% mutate(obs.freq=round(obs.freq,3), expected.freq=round(expected.freq,3))
  #
  #     DT::datatable(
  #       haplo.summaryTbl() %>%
  #         rename("Individual ID"=id), options = list(
  #           lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
  #           pageLength = 15
  #         )
  #     )
  #   })


})
