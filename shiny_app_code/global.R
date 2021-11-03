packages <- c("rstudioapi","shiny","dplyr","magrittr","tidyverse",
              "lubridate","patchwork","ggpubr",
              "ggplot2","devtools","bsplus",
              "shinydashboard","shinyWidgets",
              "slickR","svglite","hash","MMWRweek",
              "plotly","DT","shinyBS",
              "rhandsontable","foreach",
              "future","shinyFiles","shinyjs",
              "devtools","moments","zoo", "extraDistr")

##Antiquated libraries: "svglite","slickR","future"
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

##Install necessary packages from github
if(!require(lazymcmc)) devtools::install_github("jameshay218/lazymcmc", dependencies=TRUE)
library(lazymcmc)
if(!require(virosolver)) devtools::install_github("jameshay218/virosolver", dependencies=TRUE)
library(virosolver)

### GLOBAL VARIABLES ######
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Load default/example data
data(example_epi_data)
sample_epiDat <- example_epi_data

data(example_gp_partab)

simct_pars <- example_gp_partab$values
names(simct_pars) <- example_gp_partab$names

data(example_ct_data)
data(vignette_data2)

sample_ctDat <- vignette_data2

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

