library(tidyverse)
library(ggplot2)
library(lazymcmc)
library(extraDistr)
library(doParallel)

Rcpp::compileAttributes("~/Documents/GitHub/virosolver")
devtools::document("~/Documents/GitHub/virosolver")
devtools::load_all("~/Documents/GitHub/virosolver")

source("~/Documents/GitHub/ct_inference/code/priors.R")
setwd("~/Documents/viral_load_model_test/")

## MCMC management
n_clusters <- 10
cl <- makeCluster(n_clusters)
registerDoParallel(cl)

rerun <- TRUE
run_cumulative <- TRUE
save_wd <- "chains/two_samples_pos/"
run_name <- "two_samples_pos"
plot_wd <- "plots/two_samples_pos/"

## MCMC
mcmcPars <- c("iterations"=100000,"popt"=0.234,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=50000,"save_block"=1000)

## Parameters
parTab <- read.csv("~/Documents/GitHub/ct_inference/pars/parTab_test_seir.csv")
pars <- parTab$values
names(pars) <- parTab$names

## Observation times
#obs_times <- seq(50,250,by=10)
obs_times <- c(75,125,150)
ages <- 1:max(obs_times)
times <- seq(0,max(obs_times),by=1)


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
                                    convert_ct=TRUE,
                                    add_noise=rgumbel)

## Subset simulated viral loads to observation times
obs_dat <- viral_loads %>%
  filter(t %in% obs_times) %>%
  group_by(t) %>%
  sample_n(n_per_samp) %>%
  select(t, obs) %>%
  rename(ct=obs)
obs_dat <- obs_dat %>% filter(ct < 40)


res <- foreach(i=seq_along(obs_times),.packages = c("lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
  dir.create(save_wd)

  source("~/Documents/GitHub/ct_inference/code/priors.R")
  devtools::load_all("~/Documents/GitHub/virosolver")
  obs_time <- obs_times[i]
  if(run_cumulative) {
    obs_dat_tmp <- obs_dat %>% filter(t <= obs_time)
  } else {
    obs_dat_tmp <- obs_dat %>% filter(t == obs_time)
  }

  ## Epidemic cannot start after first observation time
  parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_tmp$t)
  ## Get random starting values
  startTab <- generate_viable_start_pars(parTab,obs_dat,
                                         create_posterior_func,
                                         incidence_func,
                                         prior_func_use)
  covMat <- diag(nrow(startTab))
  mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)

  output <- run_MCMC(parTab=startTab,
                      data=obs_dat_tmp,
                      INCIDENCE_FUNC=incidence_func,
                      PRIOR_FUNC = prior_func_use,
                      solve_likelihood=TRUE,
                      mcmcPars=mcmcPars,
                      filename=paste0(save_wd, run_name,"_",obs_time),
                      CREATE_POSTERIOR_FUNC=create_posterior_func,
                      mvrPars=mvrPars,
                      OPT_TUNING=0.2,
                     use_pos=TRUE)

  chain <- read.csv(output$file)
  chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]

  p_trace <- chain[,c(1,which(parTab$fixed == 0) +1)] %>%
    pivot_longer(-sampno) %>%
    ggplot() +
    geom_line(aes(x=sampno,y=value)) +
    facet_wrap(~name,scales="free_y")+
    scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),by=20000)) +
    theme_classic()

  predictions <- plot_prob_infection(chain, 100, incidence_func, 0:365,
                                     obs_dat=obs_dat_tmp,
                                     true_prob_infection = tibble(t=times,prob_infection=prob_infection))

  p1 <- predictions$plot
  model_func <- create_posterior_func(parTab,obs_dat_tmp,NULL,incidence_func,"model")
  p2 <- plot_distribution_fits(chain, obs_dat_tmp, model_func,100)
  p3 <- plot_posterior_density(chain, "R0",parTab, 0, 1000)
  p4 <- plot_posterior_density(chain, "t0",parTab, 0, 1000)

  dir.create(paste0(plot_wd,"/traces/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/predictions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/distributions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/R0_density/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/t0_density/"),recursive = TRUE)

  ggsave(paste0(plot_wd,"/traces/",run_name,"_",obs_time,"_trace.png"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/predictions/",run_name,"_",obs_time,"_predictions.png"),p1,width=7,height=4)
  ggsave(paste0(plot_wd,"distributions/",run_name,"_",obs_time,"_distributions.png"),p2,
         width=(7/5) * length(unique(obs_dat_tmp$t)),height=6)
  ggsave(paste0(plot_wd,"R0_density/",run_name,"_",obs_time,"_R0_density.png"),p3,width=5,height=3)
  ggsave(paste0(plot_wd,"t0_density/",run_name,"_",obs_time,"_t0_density.png"),p4,width=5,height=3)
}
