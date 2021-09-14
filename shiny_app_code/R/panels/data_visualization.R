## This function takes in ct value data (ct_dat), epidemic data (epi_dat),
## and selected filters (in the form of a hash table) and returns a 
## list object containing ct-specific plots in slot 2 and epi-specific 
## plots in slot 1. 
dv_plots <- function(ct_dat=sample_ctDat, epi_dat=sample_epiDat, filters=hash::hash() ) {
  ## filter data appropriately given filters
  if (!is.empty(filters)) { 
    ct_dat_long <- ctmake_long(ct_dat)
    ct_summ <- summarize_ct(ct_dat_long)
    e_dat_list <- emake_long(epi_dat)
    epi_dat_long <- e_dat_list[[1]] ## Indexing necessary
    epi_filters <- e_dat_list[[2]]
    grs_dat <- emake_grs(epi_dat_long)
    grs_plot <- plot_cases(epi_dat_long, epi_filters)
  }
  ## format data
  else { 
    ct_dat_long <- ctmake_long(ct_dat)
    ct_summ <- summarize_ct(ct_dat_long)
    e_dat_list <- emake_long(epi_dat)
    epi_dat_long <- e_dat_list[[1]] ## Indexing necessary
    epi_filters <- e_dat_list[[2]]
    grs_dat <- emake_grs(epi_dat_long)
    grs_plot <- plot_cases(epi_dat_long, epi_filters) #missing: type_filters
    comb_dat <- combine_vis_dat(ct_dat_long, ct_summ, grs_dat, epi_dat_long)
  }
  
  ct_plot_raw <- plot_ct_raw(ct_dat_long)
  ct_plot_mean <- plot_ct_mean(ct_dat_long, ct_summ)
  ct_plot_skew <- plot_ct_skew(ct_dat_long, ct_summ)
  gr_scat <- mean_gr_scatter(comb_dat) 
  pb_scat <- p_both_scatter(comb_dat)
  pm_time <- p_mean_time(comb_dat)
  ps_time <- p_skew_time(comb_dat)
  pc_conf <- p_cases_confirmed(epi_dat_long, comb_dat)
  #pgr_conf <- p_gr_confirmed(grs_dat, comb_dat) FIXME: issue w grs_dat - unresolved 
  vi_plot <- violin_plots(comb_dat, ct_dat_long)
  epi_plots <- list(pb_scat,pm_time,ps_time,pc_conf)
  ct_plots <- list(ct_plot_raw,ct_plot_mean,ct_plot_skew,vi_plot)
  plots <-list(epi_plots,ct_plots)
  return(plots)
}

## Returns data visualization main UI content
vis_content <- function(id, data) {
  ns <- NS(id) # Namespace function is necessary to module-ize server-side code
  filter_options <- get_options(data)[1]
  filter_names <- filter_req(data)
  content <- fluidRow(
     tabsetPanel(id=ns("dv_tabs"),
                 tabPanel(value="ct_panel", title="Ct View",
                   column(3,
                     wellPanel(
                            fileInput(ns("data2"),"Upload PCR-RT Ct Data (CSV)",  accept=".csv"),
                            filter_options,
                            actionButton(ns("filter_sub"),"Filter by Selected Values"),
                            downloadButton(id=ns("dp1"),"Download Plots")
                            )),
                   column(9,
                          splitLayout(cellWidths = (c("5%","90%","5%")),
                          actionButton(ns("leftSlide"),"", icon=icon("arrow-circle-left")),
                          plotOutput(ns("selectedPlot_ct")),
                          actionButton(ns("rightSlide"),"", icon=icon("arrow-circle-right"))),
                          wellPanel(textOutput(ns("plotCaptions1")))
                   )
                   ),
                 tabPanel(value="epi_panel", title="Epi View",
                   column(3,
                          wellPanel(
                                    fileInput(ns("data1"),"Upload Epidemic Data (CSV)", accept=".csv"),
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


vis_tab <- tabPanel("Data Visualization", value="vis_tab", vis_content("data_vis",india_data))


## Server-side code
load_data_vis <- function(id) {
  moduleServer(
    id,
    function(input=input, output=output, session=session) {
      rv <- reactiveValues(ct_data=sample_ctDat,
                           epi_data=sample_epiDat,
                           plots=NULL, plotlist=NULL,
                           ct_loaded=FALSE,
                           epi_loaded=FALSE)
      plot.info <-  reactiveValues(slideno=0 ,ggplot=NULL)
      #mcmc_gaus <- reactiveValues(pars=example_seir_partab)
      
      userData <- reactive({
        list(input$data1, input$data2)
      })
      
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
          rv$plotlist <- dv_plots(ct_dat=rv$ct_data, epi_dat=rv$epi_data) #FIXME: this should be dynamic
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
        browser()
        plot.info$slideno = plot.info$slideno - 1 
        indexNo <- (plot.info$slideno %% length(rv$plots)) + 1
        plot.info$ggplot = rv$plots[[indexNo]]
      })
      
      ## Dynamic Filtering 
      observeEvent(input$filter_sub, {
        epi_data <- epi_data()
        object_names <- get_options(epi_data)[2]
        filter_dict <- hash::hash()
        for(obj in object_names[[1]]) {
          filter_dict[[obj]] <- input[[obj]]
        }
        rv$plots <- dv_plots(filter_dict)
      })
  })}