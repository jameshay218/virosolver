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
        
        # Input: Checkbox if file has header ----
        #checkboxInput("header", "Header", TRUE),
        
        # Input: Select separator ----
        #radioButtons("sep", "Separator",
         #            choices = c(Comma = ",",
        #                         Semicolon = ";",
        #                         Tab = "\t"),
        #             selected = ","),
        
        # Input: Select quotes ----
        #radioButtons("quote", "Quote",
        #             choices = c(None = "",
        #                         "Double Quote" = '"',
        #                         "Single Quote" = "'"),
        #             selected = '"'),
        
        # Horizontal line ----
        #tags$hr(),
        
        # Input: Select number of rows to display ----
        #radioButtons("disp", "Display",
        #             choices = c(Head = "head",
        #                         All = "all"),
        #             selected = "head")
        
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        
        # Output: Data file ----
        tableOutput("uploaded_chains")
        
      )
    )
  )
}

runmod_tab <- tabPanel("Upload Markov Chains", mcmcup_content())