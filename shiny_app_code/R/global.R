packages <- c("shiny","dplyr","magrittr","tidyverse",
              "lubridate","patchwork","ggpubr",
              "ggplot2","devtools","bsplus",
              "shinydashboard","shinyWidgets",
              "slickR","svglite","hash","MMWRweek",
              "plotly","virosolver","DT","shinyBS",
              "rhandsontable","lazymcmc","foreach",
              "future","shinyFiles","shinyjs")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rootdir <- dirname(rstudioapi::getSourceEditorContext()$path)

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
##load packages

### GLOBAL VARIABLES ######
## FIXME: Add these to R package structure
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
sample_epiDat <- read_csv("../data/india_data_districts.csv")
sample_ctDat <- read_csv("../data/RawData_COVID_PCR Analysis_oct2020.csv") 

## load example data for testing 
data(example_gp_partab)

simct_pars <- example_gp_partab$values
names(simct_pars) <- example_gp_partab$names

data(example_ct_data)

## MCMC Parameter Names
seir_parnames <- c("viral_peak","wane_rate2",
                   "t_switch","level_switch",
                   "prob_detect","incubation",
                   "infectious")

## Generate simulated observations for testing 
sim_pars <- example_gp_partab$values
names(sim_pars) <- example_gp_partab$names
sim_test_ages <- seq(1,50,by=1)
sim_cts <- simulate_viral_loads_example(sim_test_ages, sim_pars,N=200)