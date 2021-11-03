## This function takes in ct value data (ct_dat), epidemic data (epi_dat),
## and selected filters (in the form of a hash table) and returns a 
## list object containing ct-specific plots in slot 2 and epi-specific 
## plots in slot 1. 

source(paste0(rootdir,"/helpers/plots.R"))
source(paste0(rootdir,"/global.R"))

## This function takes in either the default or user-uploaded 
## Ct value dataframe, epidemic demographic dataframe, 
## and the selected filters to be applied to each dataframe.
##
## This function creates all necessary graphs for the data 
## visualization page and returns them in list format. 
## FIXME: these plots take roughly a minute or two to load, 
## it's very slow! 
dv_plots <- function(ct_dat=sample_ctDat, epi_dat=sample_epiDat,
                     ct_filters=hash::hash(), epi_filters=hash::hash() ) {
  ## Filter data appropriately given user input 
  if (!is.empty(ct_filters)) {
    ct_dat_long <- ctmake_long(ct_dat=ct_dat, filters=ct_filters)
    browser()
    ct_summ <- summarize_ct(ct_dat_long)
    e_dat_list <- emake_long(epi_dat)
    browser()
    epi_dat_long <- e_dat_list[[1]] ## Indexing necessary
    epi_filters <- e_dat_list[[2]] 
    grs_dat <- emake_grs(epi_dat_long)
    ##FIXME: Epi tab filtering not currently in place. 
    grs_plot <- plot_cases(epi_dat_long, epi_filters)
    comb_dat <- combine_vis_dat(ct_dat_long, ct_summ, grs_dat, epi_dat_long)
  }
  ## Prepare data without filters. 
  else { 
    ct_dat_long <- ctmake_long(ct_dat)
    ct_summ <- summarize_ct(ct_dat_long)
    e_dat_list <- emake_long(epi_dat)
    epi_dat_long <- e_dat_list[[1]] ## Indexing necessary
    epi_filters <- e_dat_list[[2]]
    grs_dat <- emake_grs(epi_dat_long)
    grs_plot <- plot_cases(epi_dat_long, epi_filters)
    comb_dat <- combine_vis_dat(ct_dat_long, ct_summ, grs_dat, epi_dat_long)
  }
  
  ct_plot_raw <- plot_ct_raw(ct_dat_long)
  ct_plot_mean <- plot_ct_mean(ct_dat_long, ct_summ)
  ct_plot_skew <- plot_ct_skew(ct_dat_long, ct_summ)
  gr_scat <- mean_gr_scatter(comb_dat) 
  pb_scat <- p_both_scatter(comb_dat)
  pm_time <- p_mean_time(comb_dat)
  ps_time <- p_skew_time(comb_dat)
  #pc_conf is a list of plots (1 per gene/assay)
  pc_conf <- p_cases_confirmed(epi_dat_long, comb_dat)
  #pgr_conf <- p_gr_confirmed(grs_dat, comb_dat) FIXME: issue w grs_dat - unresolved 
  #vi_plot is a list of plots (1 per gene/assay)
  vi_plot <- violin_plots(comb_dat, ct_dat_long)
  epi_plots <- list(pb_scat,pm_time,pc_conf)
  for (ps_plot in ps_time) { epi_plots[[length(epi_plots) + 1]] <- ps_plot }
  ct_plots <- list(ct_plot_raw,ct_plot_mean,ct_plot_skew)
  for (v_plot in vi_plot) { ct_plots[[length(ct_plots) + 1]] <- v_plot }
  plots <-list(epi_plots,ct_plots)
  return(plots)
}
## This function creates the initial dropdown selector UI for 
## the user to choose which columns they would like to filter on 
## for both the Ct and Epi tab of the data visualization page. 
show_filtercandidates <- function(data,epi=FALSE, ns) {
  candidates <- colnames(data[, sapply(data, class) %in% c('character', 'factor','numeric','Date')])
  if (epi==FALSE) {id="ctfilter_candidates"}
  else {id="epifilter_candidates"}
  
  pickerInput(
    inputId = ns(id), 
    label = "Select columns for filtering:", 
    choices = sort(candidates), options = list(`actions-box` = TRUE), 
    multiple = TRUE
  )
}

## This function  creates the UI that
## allows the user to indicate which 
## columns contain Ct values (multiple selections allowed). 
## It is necessary for accurate data manipulation. 
show_ctcandidates <- function(data,ns) {
  candidates <- colnames(data[, sapply(data, class) %in% c('numeric')])
  id="assay_cols"
  pickerInput(
    inputId = ns(id), 
    label = "Please indicate which columns contain CT values:", 
    choices = sort(candidates), options = list(`actions-box` = TRUE), 
    multiple = TRUE
  )
}


## This function creates UI selector objects for 
## selected filter-on columns 
show_filtervals <- function(data, cols, epi=FALSE, ns) {
  get_options(data,cols,ns=ns)[1]
}

## This function returns data visualization main UI content
vis_content <- function(id, ct_data=sample_ctDat, epi_data=sample_epiDat) {
  ns <- NS(id) # Namespace function is necessary to modularize server-side code
  content <- fluidRow(
    useShinyjs(),
     tabsetPanel(id=ns("dv_tabs"),
                 tabPanel(value="ct_panel", title="Ct View",
                   column(3,
                     wellPanel(
                            fileInput(ns("data2"),"Upload PCR-RT CTData (CSV):",  accept=".csv"),
                            tags$hr(),
                            show_filtercandidates(ct_data,ns=ns),
                            tags$hr(),
                            show_ctcandidates(ct_data, ns=ns),
                            tags$hr(),
                            div(id="placeholder"),#Necessary for updating the DOM using this div as a reference
                            tags$hr(),
                            actionButton(ns("filter_sub"),"Filter Data"),
                            downloadButton(id=ns("dp1"),"Download Plots")  # FIXME: the download button is not fully functional 
                            )),
                   column(9,
                          splitLayout(cellWidths = (c("5%","90%","5%")),
                          div(style="display:inline-block", actionButton(ns("leftSlide"),"", icon=icon("arrow-circle-left"),align="center")),
                          plotOutput(ns("selectedPlot_ct")),
                          actionButton(ns("rightSlide"),"", icon=icon("arrow-circle-right"),align="center")),
                          wellPanel(textOutput(ns("plotCaptions1"))) # FIXME: captions have not been updated. 
                   )
                   ),
                 tabPanel(value="epi_panel", title="Epi View",
                   column(3,
                          wellPanel(
                                    fileInput(ns("data1"),"Upload Epidemic Data (CSV)", accept=".csv"),
                                    tags$hr(),
                                    show_filtercandidates(epi_data,epi=TRUE,ns=ns),
                                    tags$hr(),
                                    downloadButton(id=ns("dp2"), "Download Plots")
                          )),
                   column(9,
                          splitLayout(cellWidths = (c("5%","90%","5%")),
                          actionButton(ns("leftSlide2"),"", icon=icon("arrow-circle-left")),
                          plotOutput(ns("selectedPlot_epi")),
                          actionButton(ns("rightSlide2"),"", icon=icon("arrow-circle-right"))),
                          wellPanel(textOutput(ns("plotCaptions2")))
                   )
                 ))
  )
  return(content)
}

vis_tab <- tabPanel("Data Visualization", value="vis_tab", vis_content("data_vis", ct_data=sample_ctDat, epi_data=sample_epiDat))

## Server-side code
load_data_vis <- function(id) {
  moduleServer(
    id,
    function(input=input, output=output, session=session) {
      rv <- reactiveValues(ct_data=sample_ctDat,
                           epi_data=sample_epiDat,
                           plots=NULL, plotlist=NULL,
                           ct_loaded=FALSE,
                           epi_loaded=FALSE,
                           displayed_filters=c())
      
      
      plot.info <-  reactiveValues(slideno=0 ,ggplot=NULL)

      userData <- reactive({
        list(input$data1, input$data2)
      })
      
      # Update data used for application 
      # when a user uploads their own data. 
      observeEvent(userData(), {
        if(!is.null(input$data1)) {
          rv$epi_data <- read_csv(file=(input$data1)$datapath)
          rv$epi_loaded <- FALSE
        }
        if(!is.null(input$data2)) {
          rv$ct_data <- read_csv(file=(input$data2)$datapath)
          rv$ct_loaded <- FALSE
        } 
      })
      
      ## Data vis tab-dependent 
      observeEvent(input$dv_tabs,{
        if(rv$ct_loaded == FALSE && rv$epi_loaded == FALSE) {
          rv$plotlist <- dv_plots(ct_dat=rv$ct_data, epi_dat=rv$epi_data) 
        }
        if(input$dv_tabs == "ct_panel") {
            rv$plots <- rv$plotlist[[2]]
            plot.info$ggplot <- rv$plots[[1]]
            output$selectedPlot_ct <- renderPlot(plot.info$ggplot)
            rv$ct_loaded <- TRUE
        }
        if(input$dv_tabs == "epi_panel") {
          rv$plots <- rv$plotlist[[1]]
          plot.info$ggplot <- rv$plots[[1]]
          output$selectedPlot_epi <- renderPlot(plot.info$ggplot)
          rv$epi_loaded <- TRUE
        }
      })
      
      leftSlide <- reactive({
        list(input$leftSlide, input$leftSlide2)
      })
      
      observeEvent(leftSlide(), {
        plot.info$slideno = plot.info$slideno - 1 
        indexNo <- (plot.info$slideno %% length(rv$plots)) + 1
        plot.info$ggplot = rv$plots[[indexNo]]
      })
      
      rightSlide <- reactive({
        list(input$rightSlide, input$rightSlide2)
      })
      
      observeEvent(input$rightSlide, {
        plot.info$slideno = plot.info$slideno - 1 
        indexNo <- (plot.info$slideno %% length(rv$plots)) + 1
        plot.info$ggplot = rv$plots[[indexNo]]
      })
      
      ns <- session$ns
      
      ##Return plots
      get_plots <- reactive({ #changed from reactive
        rv$plots
      })
      
      ##Download button handler
      ## Not functional; FIXME: need to figure out how to 
      ## access reactive values for downloads. 
      output$dp1 <- downloadHandler(
        filename = function() {
          paste("viro-plots-", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file)
          get_plots()
          dev.off()
        }
      )
      

      ## Dynamic Filtering 
      filter_cols <- reactive({ #changed from reactive
        list(input$epifilter_candidates,input$ctfilter_candidates)
      })
      

      ## On filter column selection, create
      ## UI objects for selectors containing unique values for 
      ## those columns.  
      ##
      ## On unclick of a column, delete UI objects. 
      observeEvent(filter_cols(), {
        browser()
        candidates <- colnames(rv$ct_data[, sapply(rv$ct_data, class) %in% c('character', 'factor','numeric','Date')]) #place in reactive statement -- doesn't need to re-execute!
        selected_filters <- input$ctfilter_candidates
        old_filters <- rv$displayed_filters
        rv$displayed_filters <- selected_filters
        delete_filters <- old_filters[!(old_filters %in% selected_filters)]
        selected_filters <- selected_filters[!(selected_filters %in% old_filters)] #if filters are already displayed, do not create new UI object. 
        
        for (del_fil in delete_filters) {
          id <- ns(paste0(del_fil,"_options"))
          removeUI(selector=paste0("#",id), immediate = TRUE)
          removeUI(selector=sprintf("div#%s",id), immediate = TRUE)
        }
        
        for(col in selected_filters) {
          id <- ns(paste0(col,"_options"))
          if (col %in% candidates) {
            insertUI(
              selector="#placeholder",
              where="beforeBegin",
              ui = show_filtervals(rv$ct_data,cols=col,ns=session$ns)
            )
          }
        }
      } ,ignoreNULL=FALSE)
      
      ## On filter value submission, 
      ## filter dataframe and update plots. 
      observeEvent(input$filter_sub, {
        filter_objs <- rv$displayed_filters
        filter_dict <- hash::hash()
        for(obj in filter_objs) {
          obj <- paste0(obj,"_options")
          filter_val <- input[[obj]]
          if(!is.null(filter_val)) { filter_dict[[obj]] <- input[[obj]] }
        }
        
        if(input$dv_tabs == "ct_panel") { 
          rv$plotlist <- dv_plots(ct_dat=rv$ct_data, epi_dat=rv$epi_data,ct_filters=filter_dict)
          rv$plots <- rv$plotlist[[2]]
          plot.info$ggplot <- rv$plots[[1]]
          output$selectedPlot_ct <- renderPlot(plot.info$ggplot)
        }
        
        if(input$dv_tabs == "epi_panel") { rv$plotlist <- dv_plots(ct_dat=rv$ct_data, epi_dat=rv$epi_data, epi_filters=filter_dict)}
      })
  })}