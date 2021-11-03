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

#dv.packages <- c("virosolver")

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

if(!require(lazymcmc)) devtools::install_github("jameshay218/lazymcmc")
library(lazymcmc)
if(!require(virosolver)) devtools::install_github("jameshay218/virosolver")
library(virosolver)

### GLOBAL VARIABLES ######
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

load("C:/Users/anvia/Desktop/Dakotah_Testing/virosolver/data/example_epi_data.RData")
sample_epiDat <- sample.epi.df

#sample_epiDat <- read_csv("../data/india_data_districts.csv")

#sample_epiDat <- read_csv("../data/")

#sample_ctDat <- read_csv("../data/RawData_COVID_PCR Analysis_oct2020.csv")

#sample_ctDat <- read_csv("../data/")

## load example data for testing 
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

