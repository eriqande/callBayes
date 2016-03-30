library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  style="padding-top: 180px;z-index: -1;",
  # Application title
  #titlePanel("Haplotype Viewer"), 
  absolutePanel(
    style="z-index:100",
    top = 0, left = 0, right = 0,
    fixed = TRUE,
    div(
      style="padding: 6px; border-bottom: 1px solid #CCC; background: #CFCFCD",
      fluidRow( 
        column(2,column(8,offset=2,"HapPLOType",style="margin-top:16px;font-family:helvetica; font-size:40px"),
               column(2)),
        column(3,
               column(3, selectInput("selectLocus", label = "Locus:","ALL",selected = "ALL"),
                      style="padding-right: 0px; margin-top:10px;padding-left:0px; padding-right: 0px;"),
               column(1,actionButton("locusBack", label="<<"),style="margin-top:36px; padding-right: 0px;"),
               column(1,actionButton("locusFor", label=">>"),style="margin-top:36px; padding-right: 0px;"),
               column(3,selectInput("selectIndiv", label = "Individual:","ALL", selected = "ALL"),
                      style="padding-right: 0px; margin-top:10px",offset=1),
               column(1,actionButton("indivBack", label="<<"),style="margin-top:36px; padding-right: 0px;"),
               column(1,actionButton("indivFor", label=">>"),style="margin-top:36px; padding-right: 0px;"),
               column(1)
        ),
        column(3, 
               column(5, sliderInput("coverageMin",
                                     "Minimum Read Coverage",
                                     min = 0,
                                     max = 200,
                                     value = 1)),
               column(5,sliderInput("minAlleleRatio",
                                    "Minimum Allelic Depth Ratio",
                                    min = 0,
                                    max = 1,
                                    value = 0.2)),
               column(2, actionButton("updateFilter", label="Update"), style="margin-top:36px; padding-right: 0px;")),
        column(1, checkboxInput("topTwo", label = h5("Keep only the top 2 most common Haplotypes"), value = FALSE)),
        column(1, h5("Post Filtered Table:"),downloadButton('downloadData', 'Download'), align="center",
               tags$style(type='text/css', "#downloadData { vertical-align: bottom; height: 40px;margin-top:10px;font-size:20px}"),
               offset=2)
      ))),
  #hr(),
  #selectInput("selectLocus", label = "", 
  #            choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
  #           selected = 1),
  
  column(6,renderUI
         plotOutput("haplDensityPlot", height = 700,
                      dblclick = dblclickOpts(
                        id = "plotH_dblclick"),
                      brush = brushOpts(
                        id = "plotH_brush",
                        direction = "y",
                        resetOnNew = TRUE))),
  column(1,plotOutput("fracIndivPlot",height = 700)),
  column(1,plotOutput("numHapPlot",height = 700)),
  column(4, plotOutput("distPlot3",height = 700)),
  column(4,plotOutput("AlleleRatioByIndiv",
                      dblclick = dblclickOpts(
                        id = "plot_dblclick"),
                      brush = brushOpts(
                        id = "plot_brush",
                        direction = "y",
                        resetOnNew = TRUE))),
  column(4,plotOutput("readDepthByIndiv")),
  column(4,plotOutput("distPlot", 
                      dblclick = dblclickOpts(
                        id = "plot1_dblclick"),
                      brush = brushOpts(
                        id = "plot1_brush",
                        direction = "y",
                        resetOnNew = TRUE))),
  column(4,plotOutput("distPlot2")),
  column(4,plotOutput("histHap")),
  column(4,plotOutput("PairWiseHap")),
  column(4, DT::dataTableOutput('haploTbl'),style="border-right:2px solid grey;"),
  column(6, DT::dataTableOutput('haploSummary'),
         offset=2),
  titlePanel("", windowTitle = "HapPLOType: a view to your haplotypes")
  
))
