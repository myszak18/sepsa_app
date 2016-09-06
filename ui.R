library(shiny)

# Define UI for random distribution application 
shinyUI(fluidPage(
    
  # Application title
  titlePanel("Sepsis"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the
  # br() element to introduce extra vertical spacing
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("organ", "Organs:",choices=
                  c("liver","kidney","heart","lungs")),
      selectInput("GO_type", "GO analysis",
                  c("none", "MF","BP","CC")),
      actionButton("goButton", "Submit"),
      
      br(),
      br(),
    downloadButton('downloadData', 'Download')
    ),
    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Table", dataTableOutput("table")),
        tabPanel("Plot", plotOutput("plot") ),
        tabPanel("Plot GO up", plotOutput("plot_go_up",click = "plot_click"), dataTableOutput("info") ),
        tabPanel("Plot GO down", plotOutput("plot_go_down",click = "plot_click"), dataTableOutput("info") )
      )
    )
    )
  )
)

