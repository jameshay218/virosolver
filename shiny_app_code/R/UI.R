setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rootdir <- dirname(rstudioapi::getSourceEditorContext()$path)

## Load in helper functions
source(paste0(rootdir,"/helpers/plots.R"))
source(paste0(rootdir,"/helpers/mcmc_helpers.R"))
source(paste0(rootdir,"/helpers/vk_helpers.R"))
source(paste0(rootdir,"/helpers/file_management.R"))

## Load in UI panels
source(paste0(rootdir,"/panels/data_visualization.R"))
source(paste0(rootdir,"/panels/overview.R"))
source(paste0(rootdir,"/panels/viral_kinetics.R"))
source(paste0(rootdir,"/panels/mcmc_params.R"))
source(paste0(rootdir,"/panels/run_model.R"))
source(paste0(rootdir,"/panels/primary_output.R"))



viro_UI <- shinyUI ({
  fluidPage(
    navbarPage(id="viro_nav",
               "Virosolver",
               overview_tab,
               vis_tab, 
               vk_tab,
               mcmc_tab,
               runmod_tab,
               out_tab
    ) 
  )
})