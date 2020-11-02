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
run_cumulative <- FALSE
save_wd <- "chains/pos_exp/"
run_name <- "pos_exp"
plot_wd <- "plots/pos_exp/"

## MCMC
mcmcPars <- c("iterations"=20000,"popt"=0.234,"opt_freq"=1000,
              "thin"=10,"adaptive_period"=10000,"save_block"=1000)

## Parameters
parTab <- read.csv("~/Documents/GitHub/ct_inference/pars/parTab_test_seir.csv")
parTab[parTab$names == "overall_prob","values"] <- 1
parTab[parTab$names == "overall_prob","fixed"] <- 1
pars <- parTab$values
names(pars) <- parTab$names

## Observation times
obs_times <- seq(50,250,by=25)
ages <- 1:max(obs_times)
ages <- 1:35
times <- seq(0,max(obs_times),by=1)


## Specify process models
prior_func_use <- prior_func_hinge_exp
incidence_func_sim <- solveSEIRModel_rlsoda_wrapper
incidence_func <- exponential_growth_model

## Number sampled per time
n_per_samp <- 5000
n_overall <- length(obs_times)*n_per_samp

## Probability of infection

## Simulate probability of infection
prob_infection <- incidence_func_sim(pars, times)
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

  obs_dat_tmp1 <- obs_dat_tmp
  obs_dat_tmp$t <- max(ages)

  ## Get random starting values
  startTab <- parTab
  #startTab <- generate_viable_start_pars(parTab,obs_dat,
  #                                       create_posterior_func,
  #                                       incidence_func,
  #                                       prior_func_use)

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

  max_obs_t <- max(obs_dat_tmp1$t)
  tmp_times <- (max_obs_t - length(ages)):max_obs_t


  prob1 <- prob_infection[which(times %in% tmp_times)]
  prob2 <- prob_infection/sum(prob1)
  predictions <- plot_prob_infection(chain, 100, incidence_func, seq(max(obs_dat_tmp1)-max(ages),max(obs_dat_tmp1$t)),
                                     obs_dat=obs_dat_tmp1,
                                     true_prob_infection = tibble(t=times,prob_infection=prob2),
                                     tshift=0)

  p1 <- predictions$plot
  model_func <- create_posterior_func(parTab,obs_dat_tmp,NULL,incidence_func,"model")
  p2 <- plot_distribution_fits(chain, obs_dat_tmp, model_func,100)
  p3 <- plot_posterior_density(chain, "beta",parTab, 0, 1000)
  p4 <- plot_posterior_density(chain, "overall_prob",parTab, 0, 1000)

  dir.create(paste0(plot_wd,"/traces/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/predictions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/distributions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/exp/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/overall_prob/"),recursive = TRUE)

  ggsave(paste0(plot_wd,"/traces/",run_name,"_",obs_time,"_trace.png"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/predictions/",run_name,"_",obs_time,"_predictions.png"),p1,width=7,height=4)
  ggsave(paste0(plot_wd,"distributions/",run_name,"_",obs_time,"_distributions.png"),p2,
         width=(7/5) * length(unique(obs_dat_tmp$t)),height=6)
  ggsave(paste0(plot_wd,"exp/",run_name,"_",obs_time,"_exp.png"),p3,width=5,height=3)
  ggsave(paste0(plot_wd,"overall_prob/",run_name,"_",obs_time,"_overall_prob.png"),p4,width=5,height=3)
}
