BiocManager::install(c("flowCore", "flowViz"))

library(shiny)
library(flowCore)
library(flowViz)

ui <- fluidPage(
  plotOutput("plot", click = "plot_click"),
  verbatimTextOutput("gate_info")
)

server <- function(input, output, session) {
  # Load FCS data
  fcs <- read.FCS("path/to/your/file.fcs")
  
  # Reactive gate
  gate <- reactiveVal(NULL)
  
  output$plot <- renderPlot({
    plot(fcs, main = "Click to Gate")
    if (!is.null(gate())) {
      rect(gate()[1], gate()[2], gate()[3], gate()[4], border = "red")
    }
  })
  
  observeEvent(input$plot_click, {
    click <- input$plot_click
    if (is.null(gate())) {
      gate(c(click$x, click$y, click$x + 100, click$y + 100))  # Example gate
    }
  })
  
  output$gate_info <- renderPrint({
    gate()
  })
}

shinyApp(ui, server)

