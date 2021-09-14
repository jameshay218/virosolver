mcmcup_content <- function() {
  fluidPage(
    # App title ----
    titlePanel("Uploading Files"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Select a file ----
        fileInput("mcmc_chains_files", "Upload Markov Chains (.csv)",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        #shinyDirButton('mcmc_chain_folder', 'Select a folder', 'Please select a folder', FALSE),
        
        # Horizontal line ----
        tags$hr(),
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        
        # Output: Data file ----
        #tableOutput("mcmc_chain_trace")
        plotOutput("mcmc_chain_trace",height=1000)
        #tableView("uploaded_chains")
        
      )
    )
  )
}

runmod_tab <- tabPanel("Upload Markov Chains", mcmcup_content())