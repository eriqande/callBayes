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
                column(2, sliderInput("bins",
                                  "Read Coverage filter",
                                  min = 1,
                                  max = 500,
                                  value = 30)
                      )))),
      #hr(),
      #selectInput("selectLocus", label = "", 
      #            choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
       #           selected = 1),

    column(6,plotOutput("distPlot")),
    column(6,plotOutput("distPlot1"))
))
