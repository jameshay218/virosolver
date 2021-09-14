output_content <- function() { 
  fluidRow(
    mainPanel(
      plotOutput("dist_plot"),
      plotOutput("trace_plot")
    )
  )
}

out_tab <- tabPanel("Output",output_content())