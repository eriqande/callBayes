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
        column(2,column(8,offset=2,"HapPLOType",style="margin-top:15%;font-family:helvetica; font-size:200%"),
               column(2)),
        column(3,
               column(12, 
                      column(1,"Locus:",style="margin-top:10%;"),
                      column(6,selectInput("selectLocus", label ="","ALL",selected = "ALL"),
                             style="padding-right: 0px; margin-top:0.3%;padding-left:10%; padding-right: 0px;"),
                      column(2,actionButton("locusBack", label="<"), 
                             style="margin-top:6%; padding-right:0%"),
                      column(2,actionButton("locusFor", label=">"),
                             style="padding-right: 0px; margin-top:6%;padding-left:0%; margin-left: 0px; padding-right: 0px;"),
                      style="padding-right: 0px;padding-left:0px; padding-right: 0px;margin: -3% 0 0 0;"),
               #column(1,actionButton("locusBack", label="<<"),style="margin-top:10%; padding-right: 0px;"),
               #column(1,actionButton("locusFor", label=">>"),style="margin-top:10%; padding-right: 0px;"),
               column(12, 
                      column(1,"Indiv:",style="margin-top:10%;"),
                      column(6,selectInput("selectIndiv", label = "","ALL", selected = "ALL"),
                             style="padding-right: 0px; margin-top:0.3%;padding-left:10%; padding-right: 0px;"),
                      column(2,actionButton("indivBack", label="<"),
                             style="margin-top:6%; padding-right: 0px;"),
                      column(2,actionButton("indivFor", label=">"),
                             style="padding-right: 0px; margin-top:6%;padding-left:0%; margin-left: 0px; padding-right: 0px;"),
                      style="padding: 0 0 0 0; margin: -8% 0 0 0;"
                      ) 
        ),
        column(3,
               column(6,"Min read coverage:"),
               column(6,"Min allelic ratio:"),
               column(6, sliderInput("coverageMin",
                                     "",
                                     min = 0,
                                     max = 200,
                                     value = 1)),
               column(6,sliderInput("minAlleleRatio",
                                    "",
                                    min = 0,
                                    max = 1,
                                    value = 0.2))),
        column(1, checkboxInput("topTwo", label = "Keep only the top 2 most common Haplotypes", value = FALSE)),
        column(1, actionButton("updateFilter", label="Update"), style="margin-top:3%; padding-right: 0px;"
        ),
        column(2, h5("Post Filtered Table:"),downloadButton('downloadData', 'Download'), align="center",
               tags$style(type='text/css', "#downloadData { vertical-align: bottom; height: 40px;margin-top:10px;font-size:20px}"))
      ))),
  #hr(),
  #selectInput("selectLocus", label = "", 
  #            choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
  #           selected = 1),
  #column(12, h1("")),
  fluidRow( 
  column(5,plotOutput("haplDensityPlot",height="auto",
                      dblclick = dblclickOpts(
                        id = "plotH_dblclick"),
                      brush = brushOpts(
                        id = "plotH_brush",
                        direction = "y",
                        resetOnNew = TRUE))),
  column(2,plotOutput("numHapPlot",height="auto")),
  column(2,plotOutput("fracIndivPlot",height="auto")),
  column(3, plotOutput("readDepthPerLocus",height="auto"))),
  
  column(12, h1("")),
  br(),
  fluidRow( 
  column(5,plotOutput("AlleleRatioByIndiv",height="auto",
                      dblclick = dblclickOpts(
                        id = "plot_dblclick"),
                      brush = brushOpts(
                        id = "plot_brush",
                        direction = "y",
                        resetOnNew = TRUE))),
  column(2,plotOutput("fracHaploPlot",height="auto")),
  column(2,plotOutput("meanReadDepthByIndiv",height="auto")),
  column(3,plotOutput("readDepthByIndiv",height="auto"))),
#   column(4,plotOutput("distPlot", 
#                       dblclick = dblclickOpts(
#                         id = "plot1_dblclick"),
#                       brush = brushOpts(
#                         id = "plot1_brush",
#                         direction = "y",
#                         resetOnNew = TRUE))),
  column(4,plotOutput("hapSeq",height="auto")),
  column(4,plotOutput("histHap",height="auto")),
  column(4,plotOutput("PairWiseHap",height="auto")),
fluidRow(
  div(
    style="padding: 10px; border-bottom: 8px solid white; background: white"
  ),
  column(5, DT::dataTableOutput('haploTbl')),#,style="border-right:2px solid grey;"),
  column(5, DT::dataTableOutput('haploFreqTbl'), offset=2),
  column(12),
  div(
    style="padding: 10px; border-bottom: 8px solid white; background: white"
  ),
  column(4, DT::dataTableOutput('haploSummary'))),
  titlePanel("", windowTitle = "HapPLOType: a view to your haplotypes")
  
))
