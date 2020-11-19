library(tidyverse)
library(ggplot2)
library(lazymcmc)
library(extraDistr)
Rcpp::compileAttributes("~/Documents/GitHub/virosolver")
devtools::document("~/Documents/GitHub/virosolver")
devtools::load_all("~/Documents/GitHub/virosolver")

source("~/Documents/GitHub/ct_inference/code/priors.R")
setwd("~/Documents/viral_load_model_test/")

rerun <- TRUE

## Parameters
parTab <- read.csv("~/Documents/GitHub/ct_inference/pars/parTab_test_seir.csv")
pars <- parTab$values
names(pars) <- parTab$names

## Observation times
obs_times <- c(120)
ages <- 1:max(obs_times)
times <- seq(0,max(obs_times),by=1)

## Epidemic cannot start after first observation time
parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_times)

## Specify process models
prior_func_use <- prior_func_hinge2
incidence_func <- solveSEIRModel_rlsoda_wrapper

## Number sampled per time
n_per_samp <- 2000
n_overall <- length(obs_times)*n_per_samp

## Probability of infection
pars["sigma"] <- 1/pars["incubation"]
pars["gamma"] <- 1/pars["infectious"]
pars["beta"] <- pars["R0"]*pars["gamma"]

## Simulate probability of infection
prob_infection <- incidence_func(pars, times)
## Simulate infection times
inf_times <- simulate_infection_times(n_per_samp,prob_infection)

## Simulate viral loads/ct values
viral_loads <- simulate_viral_loads(inf_times, times, pars,
                                    additional_detect_process = TRUE,
                                    convert_vl=FALSE,
                                    add_noise=rgumbel)

## Subset simulated viral loads to observation times
obs_dat <- viral_loads %>%
  filter(t %in% obs_times) %>%
  group_by(t) %>%
  sample_n(n_per_samp) %>%
  select(t, obs) %>%
  rename(ct=obs)

## MCMC
mcmcPars <- c("iterations"=100000,"popt"=0.234,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=50000,"save_block"=1000)

## Get random starting values
startTab <- generate_viable_start_pars(parTab,obs_dat,
                                       create_posterior_func,
                                       incidence_func,
                                       prior_func_use)
covMat <- diag(nrow(startTab))
mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)

## Run multivariate proposal MCMC
if(rerun){
  output <- run_MCMC(parTab=startTab,
                      data=obs_dat,
                      INCIDENCE_FUNC=incidence_func,
                      PRIOR_FUNC = prior_func_use,
                      solve_likelihood=TRUE,
                      mcmcPars=mcmcPars,
                      filename="test",
                      CREATE_POSTERIOR_FUNC=create_posterior_func,
                      mvrPars=mvrPars,
                      OPT_TUNING=0.2)

  chain <- read.csv(output$file)
} else {
  chain <- read.csv("test_multivariate_chain.csv")
}
chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]

p_trace <- chain[,c(1,which(parTab$fixed == 0) +1)] %>%
  pivot_longer(-sampno) %>%
  ggplot() +
  geom_line(aes(x=sampno,y=value)) +
  facet_wrap(~name,scales="free_y")+
  scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),by=20000)) +
  theme_classic()

predictions <- plot_prob_infection(chain, 100, incidence_func, 0:365,
                                   obs_dat=obs_dat,
                                   true_prob_infection = tibble(t=times,prob_infection=prob_infection))
p1 <- predictions$plot

model_func <- create_posterior_func(parTab,obs_dat,NULL,incidence_func,"model")
p2 <- plot_distribution_fits(chain, obs_dat, model_func,100)
p3 <- plot_posterior_density(chain, "R0",parTab, 0, 1000)
