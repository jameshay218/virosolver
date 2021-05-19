########################################
## SCRIPT TO SIMULATE LINE LIST DATA AND CT VALUES FROM AN SEIR PROCESS
## James Hay <jameshay218@gmail.com>
## 18th May 2021
## Note that this script requires the repo at `virosolver_paper` to be downloaded, as this contains auxiliary functions
########################################
## 1. Headers
########################################
## Install these packages
library(tidyverse)
library(ggplot2)
library(extraDistr)
library(patchwork)
library(ggthemes)
library(odin) ## install from CRAN
library(doParallel)
library(fitdistrplus)
library(virosolver) ## install from devtools::install_github("jameshay218/virosolver")

set.seed(1)

HOME_WD <- "~/Documents/GitHub/"

## Where to perform the simulations
MAIN_WD <- paste0(HOME_WD,"/virosolver/inst/extdata")
setwd(MAIN_WD)

########################################
## 2. Model parameters and simulation settings
########################################
## Model parameters
#example_seir_partab <- read.csv("~/Documents/GitHub/virosolver_paper/pars/massachusetts/partab_seir_model.csv")
#example_seir_partab[example_seir_partab$names == "I0","values"] <- 1/10000
#example_seir_partab[example_seir_partab$names == "R0","values"] <- 1.8
#save(example_seir_partab,file="~/Documents/GitHub/virosolver/data/example_seir_partab.RData")
data(example_seir_partab)
pars <- example_seir_partab$values
names(pars) <- example_seir_partab$names

## IMPORTANT CHECK
## - Check that the assumed viral kinetics are in line
##   with your data. This means things like peak Ct value,
##   waning rate, level of the "plateau" phase, and the 
##   limit of detection
test_ages <- seq(0,50,by=1)
cts <- viral_load_func(pars, test_ages)
prop_detect <- prop_detectable(test_ages[test_ages > 0],pars, cts[test_ages > 0])
p1 <- ggplot(data.frame(ct=cts,t=test_ages)) + 
  geom_line(aes(x=t,y=ct)) + 
  scale_y_continuous(trans="reverse",
                     limits=c(40,0)) +
  ylab("Modal Ct value") +
  xlab("Days since infection")
p2 <- ggplot(data.frame(p=c(0,prop_detect),t=test_ages)) + 
  geom_line(aes(x=t,y=p)) + 
  ylab("Proportion of infections still detectable") +
  xlab("Days since infection")
p1/p2

## Simulation parameters
population_n <- 100000
## Over the course of 250 days
times <- 0:250

## Sampling procedure - how often do we take samples from the population?
sampling_frequency <- 14
## How many samples do we take on each sample day?
sampling_number <- 1000


########################################
## 3. Full simulated line list
########################################
## Deterministic SEIR model
epidemic_process <- simulate_seir_process(pars,times,N=population_n)
## Widen solution and extract key variables
res <- epidemic_process$seir_outputs %>% pivot_wider(names_from="variable",values_from="value")
res <- res %>% rename(step=time)
incidence <- epidemic_process$per_cap_incidence
overall_prob <- epidemic_process$overall_prob_infection

## Basic plot, showing the proportion of the population infected each day.
## This is the epidemic incidence curve that you're trying to estimate.
## So you can substitute this object for **any** incidence curve that
## you want to estimate.
plot(incidence)

## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
## These entries give the **true** infection (or lack of infection) timings for each individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(incidence*population_n,times=times,population_n=population_n)

########################################
## 4. Simulate observation process
########################################
## This bit of code generates linelist data for some specified observation process.
## Here, we simulate sampling 2000 people at random from the population every 14 days
## arguments changed above
sample_probs <- c(rep(0, sampling_frequency-1),sampling_number/population_n)
sample_probs <- rep(sample_probs, length(times)/sampling_frequency +1)
sample_probs <- sample_probs[1:length(times)]
frac_report <- tibble(t=times,prob=sample_probs)
frac_report <- frac_report %>% filter(t >= 50 & t <= 150)

## frac_report is a table, giving the proportion (prob) of the population
## sampled on day t
head(frac_report)

## This function takes the complete linelist and sub-samples a proportion
## of the entire population on the days specified by `frac_report`
observed_linelist <- simulate_reporting(complete_linelist, 
                                        frac_report=NULL,
                                        timevarying_prob=frac_report,
                                        solve_times=times, 
                                        symptomatic=FALSE)

## Simulate viral load/Ct value for each person in the line list
simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals,
                                                      kinetics_pars=pars)

## Clean data to form expected by virosolver
obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
  rename(t = sampled_time, ct=ct_obs) %>% arrange(t)

example_ct_data <- obs_dat
#save(example_ct_data,file="~/Documents/GitHub/virosolver/data/example_ct_data.RData")

## Save SEIR plots
## Ct distribution plot
p_dat <- ggplot(obs_dat %>% filter(ct < pars["intercept"])) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(aes(x=t,y=ct),size=0.1,width=2,height=0) + 
  scale_y_continuous(trans="reverse") +
  theme_bw() +
  scale_x_continuous(limits=c(min(times),max(times)+50)) +
  ylab("Ct value") +
  xlab("Observation time")
p_dat

example_seir_incidence <- tibble(t=0:250,prob_infection=incidence)
save(example_seir_incidence,file="~/Documents/GitHub/virosolver/data/example_seir_incidence.RData")
