viro_server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) #Expand max input size
  
  ## Initializes reactive values 
  vk <- reactiveValues(sliderValues=example_gp_partab, vk_graph=NULL, vk_cts=sim_cts)
  mcmc_1xs <- reactiveValues(seir_pars=example_seir_partab)
  mcmc_mxs <- reactiveValues(gp_pars=example_gp_partab[1:14,])

  load_data_vis("data_vis")
  
  ## Viral Kinetics Page Elements
  listen_sliders <- reactive({ list(
    input$overall_prob,
    input$tshift,
    input$desired_mode,
    input$viral_peak,
    input$t_switch,
    input$obs_sd,
    input$sd_mod,
    input$sd_mod_wane,
    input$intercept,
    input$t_switch,
    input$level_switch,
    input$wane_rate2,
    input$prob_detect
  )})
  
  output$vk_iplot <- renderPlotly(plot_vk(vk$sliderValues, vk$vk_cts))
  output$test_sliders <- renderText({vk$sliderValues$values[vk$sliderValues$names=="overall_prob"]})
  
  observeEvent(listen_sliders(),{#FIXME: move to function SWITCHED FROM observe()
    vk$sliderValues$values[vk$sliderValues$names=="overall_prob"] <- input$overall_prob
    vk$sliderValues$values[vk$sliderValues$names=="tshift"] <- input$tshift
    vk$sliderValues$values[vk$sliderValues$names=="desired_mode"] <- input$desired_mode
    vk$sliderValues$values[vk$sliderValues$names=="viral_peak"] <- input$viral_peak
    vk$sliderValues$values[vk$sliderValues$names=="t_switch"] <- input$t_switch
    vk$sliderValues$values[vk$sliderValues$names=="obs_sd"] <- input$obs_sd
    vk$sliderValues$values[vk$sliderValues$names=="sd_mod"] <- input$sd_mod
    vk$sliderValues$values[vk$sliderValues$names=="sd_mod_wane"] <- input$sd_mod_wane
    vk$sliderValues$values[vk$sliderValues$names=="intercept"] <- input$intercept
    vk$sliderValues$values[vk$sliderValues$names=="t_switch"] <- input$t_switch
    vk$sliderValues$values[vk$sliderValues$names=="level_switch"] <- input$level_switch
    vk$sliderValues$values[vk$sliderValues$names=="wane_rate2"] <- input$wane_rate2
    vk$sliderValues$values[vk$sliderValues$names=="prob_detect"] <- input$prob_detect
  })
  
  ##Viral kinetics uploads must currently include "age" and "ct" column explicitly
  observeEvent(input$vk_data, {
    vk$vk_cts <- read_csv(file=(input$vk_data)$datapath)
    ## NOTE: plotly allows for better graphic interaction and downloading, but 
    ## causes EXTREMELY slow loading with many datapoints. 
    #output$vk_iplot <- renderPlotly(plot_vk(vk$sliderValues, vk$vk_cts))
  })
  
  
  ## MCMC Parameters 
  
  ## Single Cross Section
  ## Render editable, readable table for user SEIR 
  ## parameter input 
  output$seir_pars <- renderRHandsontable({
    rhandsontable(mcmc_1xs$seir_pars)
  })
  
  observeEvent(input$singlexs_sub, {
    seir_partab <- hot_to_r(input$seir_pars) #converts table user input to dataframe
    withProgress(message="Running MCMC Framework", value=0, { # FIXME: Progress bar NOT functional 
      mcmc_info <- single_xs(time=111,ct_data=example_ct_data, #FIXME: time is user-input, needs to be added to UI!
                             seir_partab=seir_partab)
    })
    output$dist_plot <- mcmc_info[[1]]
    output$trace_plot <- mcmc_info[[2]]
  })
  
  ## Multiple Cross Sections
  output$gp_pars <- renderRHandsontable({
    rhandsontable(mcmc_mxs$gp_pars)
  })
  
}