
textInputRow<-function (inputId, label, value = "") 
{
  chunk <- div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small")
   )
  chunk <- c(
    div(style="display:inline-block",
        tags$label(label, `for` = inputId), 
        tags$input(id = paste0(inputId,"_f"), type="checkbox")
    ), chunk)
  chunk
}

## FIXME: Automate this, this is so ugly. 
## FIXME: Implement hash or named list for variable descriptions AND values,
## same for data vis! 
value_view <- function() {
  opts <- c()
  for (par in seir_parnames) {
    opt <- fluidRow(
        textInput(paste(par), paste(par),value=100000), #temp values
        checkboxInput(paste0(par,"f"), label=NULL, value=FALSE)
    )
    opts <- c(opts, opt)
  }
  content <- fluidPage(
    column(8, h2("Parameters for MCMC Calc.")),
    column(4, h2("Fixed")),
    opts,
    actionButton("singlexs_sub","Compute MCMC Chains")
    )
}

table_view <- function() {
  content <- fluidRow (
    h2("Parameters for MCMC Chain calculation"),
    rHandsontableOutput("seir_pars"),
    actionButton("singlexs_sub","Compute MCMC Chains")
  )
}

single_cross <- function(calculating=FALSE, table_view=TRUE) {
  if (table_view) table_view()
  else value_view()
}

table_view_multi <- function() { #FIXME: combine functions, this is unnecessary
  content <- fluidRow (
    h2("Parameters for MCMC Chain calculation"),
    rHandsontableOutput("gp_pars"),
    actionButton("multixs_sub","Compute MCMC Chains")
  )
}

multi_cross <- function(calculating=FALSE, table_view=TRUE) {
  table_view_multi()
}

mcmc_content <- function(table=TRUE) {
  content <- fluidRow(
  wellPanel(
    tabsetPanel(type="tabs",
                tabPanel("Single Cross Section",
                         tabsetPanel(
                           tabPanel("Upload SEIR Params",fileInput("seir_pars_USER",
                                                                   "Upload csv containing SEIR parameters", 
                                                                   accept=".csv")),
                           tabPanel("Table View", single_cross()),
                           tabPanel("Value View", single_cross(table_view=FALSE))
                         )),
                tabPanel("Multiple Cross Sections",
                         tabsetPanel(
                           tabPanel("Upload GP Params",fileInput("gp_pars_USER",
                                                                   "Upload csv containing GP parameters", 
                                                                   accept=".csv")),
                           tabPanel("Table View", multi_cross()))#,
                  )
                #tabPanel("Test multi",tabsetPanel(tabPanel("take1"),tabPanel("take2")))
    )
  )
  )
  #content <- 
  #  wellPanel(
  #    textOutput("Calculating MCMC Chains...")
  #  )
}

## generate prior function for SEIR input
prior_func_seir <- function(partable=FALSE) {
  
}

mcmc_tab <- tabPanel("MCMC", mcmc_content())