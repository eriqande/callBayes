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
            column(3,titlePanel("  Haplotype Viewer")),
                column(2,uiOutput("uiLocus")),
               column(2, uiOutput("uiIndiv")),
                column(2, sliderInput("coverageMin",
                                  "Read Coverage filter",
                                  min = 0,
                                  max = 200,
                                  value = 1)
                      )))),
      #hr(),
      #selectInput("selectLocus", label = "", 
      #            choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
       #           selected = 1),

    plotOutput("haplDensityPlot", height = 700,
               dblclick = dblclickOpts(
                 id = "plotH_dblclick"),
               brush = brushOpts(
                 id = "plotH_brush",
                 direction = "y",
                 resetOnNew = TRUE)),
    column(6,plotOutput("distPlot",
                        dblclick = dblclickOpts(
                          id = "plot_dblclick"),
                        brush = brushOpts(
                          id = "plot_brush",
                          direction = "y",
                          resetOnNew = TRUE))),
    column(6,plotOutput("distPlot1", 
                        dblclick = dblclickOpts(
                          id = "plot1_dblclick"),
                        brush = brushOpts(
                          id = "plot1_brush",
                          direction = "y",
                          resetOnNew = TRUE))),
    column(4,plotOutput("distPlot2")),
  column(4,plotOutput("distPlot3")),
  column(4,plotOutput("distPlot4")),
  DT::dataTableOutput('haploTbl')
    
))
