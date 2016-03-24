library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  style="padding-top: 180px;z-index: -1;",
  # Application title
  #titlePanel("Haplotype Viewer"),
  
  # Sidebar with a slider input for the number of bins
    absolutePanel(
      style="z-index:100",
      top = 0, left = 0, right = 0,
      fixed = TRUE,
        div(
          style="padding: 8px; border-bottom: 1px solid #CCC; background: #CFCFCD",
          fluidRow( 
            column(2,titlePanel("Haplotype Viewer"), offset = 1),
                column(1,uiOutput("uiLocus")),
               column(1, uiOutput("uiIndiv")),
                column(2, sliderInput("coverageMin",
                                  "Read Coverage Filter",
                                  min = 0,
                                  max = 200,
                                  value = 1)),
                       column(2, sliderInput("minAlleleRatio",
                                             "Minimum Allelic Depth Ratio (to hap1)",
                                             min = 0,
                                             max = 1,
                                             value = 0.5)),
            column(1, checkboxInput("topTwo", label = h5("Keep only top 2 Haplotypes"), value = FALSE)),
            column(1, h5("Post Filtered Table:"),downloadButton('downloadData', 'Download'), align="center",
                   tags$style(type='text/css', "#downloadData { vertical-align: bottom; height: 40px;margin-top:10px;font-size:20px}"))
                      ))),
      #hr(),
      #selectInput("selectLocus", label = "", 
      #            choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
       #           selected = 1),

    column(8,plotOutput("haplDensityPlot", height = 700,
               dblclick = dblclickOpts(
                 id = "plotH_dblclick"),
               brush = brushOpts(
                 id = "plotH_brush",
                 direction = "y",
                 resetOnNew = TRUE))),
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
  column(4, DT::dataTableOutput('haploTbl')),
  column(6, DT::dataTableOutput('haploSummary'),
         offset=2)
    
))
