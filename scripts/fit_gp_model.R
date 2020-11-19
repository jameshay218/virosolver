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
parTab[parTab$names == "beta","fixed"] <- 1
parTab[parTab$names == "rho","fixed"] <- 1
parTab[parTab$names == "rho","values"] <- 0.05
parTab[parTab$names == "nu","fixed"] <- 0
parTab[parTab$names == "nu","values"] <- 3
parTab[parTab$names == "t0","values"] <- 0
parTab[parTab$names == "t0","fixed"] <- 1
pars <- parTab$values
names(pars) <- parTab$names

## Observation times
obs_times <- seq(25,100,by=25)
obs_times <- 75
#obs_times <- c(35, 70,105)

ages <- 1:max(obs_times)
times <- seq(0,max(obs_times),by=1)
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))

## Epidemic cannot start after first observation time
parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_times)

## Specify process models
prior_func_use <- prior_func_hinge_gp
incidence_func <- solveSEIRModel_rlsoda_wrapper

## Number sampled per time
n_per_samp <- 5000
n_overall <- length(obs_times)*n_per_samp

## Probability of infection
pars["sigma"] <- 1/pars["incubation"]
pars["gamma"] <- 1/pars["infectious"]
pars["beta"] <- pars["R0"]*pars["gamma"]

## Simulate probability of infection
prob_infection <- incidence_func(pars, times)

#mat <- matrix(rep(times, each=length(times)),ncol=length(times))
#dist <- abs(apply(mat, 2, function(x) x-times))
#mus <- rep(0, length(times))
#nusq <- 10
#rhosq <- 0.001
#K <- nusq * exp(-rhosq * (dist)^2)
#diag(K) <- diag(K) + 0.01
#k <- mvtnorm::rmvnorm(1, mean=mus, sigma=K)[1,]
#ps <- 1/(1+exp(-k))
#prob_infection <- (ps/sum(ps))*pars["overall_prob"]
#plot(prob_infection)
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

## Set up GP parameters
tmp <- prob_infection[(0:max(obs_dat$t))+1]
tmp <- pmax(0.0000001, tmp)
true_probs <- log(tmp/(1-tmp))
names <- "prob"
fixed <- 0
lower_bound <- 0
upper_bound <- 1
steps <- 0.1
lower_start <- 0
upper_start <- 0.0001

parTab_probs <- data.frame(values=true_probs,names=names,fixed=0,
                           lower_bound=-20,upper_bound=20,steps=0.1,
                           lower_start=-2,upper_start=2)
parTab <- rbind(parTab, parTab_probs)

pars <- parTab$values
names(pars) <- parTab$names


## MCMC
mcmcPars1 <- c("iterations"=50000,"popt"=0.44,"opt_freq"=5000,
               "thin"=10,"adaptive_period"=20000,"save_block"=1000)
mcmcPars2 <- c("iterations"=100000,"popt"=0.234,"opt_freq"=5000,
              "thin"=10,"adaptive_period"=20000,"save_block"=1000)

## Get random starting values
startTab <- generate_viable_start_pars(parTab,obs_dat,
                                       create_posterior_func,
                                       gaussian_process_model,
                                       prior_func_use)
covMat <- diag(nrow(startTab))
mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)

f <- create_posterior_func(startTab, data=obs_dat, PRIOR_FUNC=prior_func_hinge_gp,
                           INCIDENCE_FUNC=gaussian_process_model,t_dist=t_dist)
f(pars)
## Run multivariate proposal MCMC
if(rerun){
  output <- run_MCMC(parTab=startTab,
                      data=obs_dat,
                      INCIDENCE_FUNC=gaussian_process_model,
                      PRIOR_FUNC = prior_func_hinge_gp,
                      solve_likelihood=TRUE,
                      mcmcPars=mcmcPars1,
                      filename="test1",
                      CREATE_POSTERIOR_FUNC=create_posterior_func,
                      mvrPars=NULL,
                      OPT_TUNING=0.2,
                     t_dist=t_dist)
  chain <- read.csv(output$file)

  bestPars <- get_best_pars(chain)
  chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
  covMat <- cov(chain)

  covMat <- diag(nrow(startTab))
  mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)

  startTab$values <- bestPars

  output1 <- run_MCMC(parTab=startTab,
                     data=obs_dat,
                     INCIDENCE_FUNC=gaussian_process_model,
                     PRIOR_FUNC = prior_func_hinge_gp,
                     solve_likelihood=TRUE,
                     mcmcPars=mcmcPars2,
                     filename="test1",
                     CREATE_POSTERIOR_FUNC=create_posterior_func,
                     mvrPars=mvrPars,
                     OPT_TUNING=0.2,
                     t_dist=t_dist)
  chain <- read.csv(output1$file)
} else {
  chain <- read.csv("test_multivariate_chain.csv")
}
chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]

samps <- sample(unique(chain$sampno), 10)
tmp_pars <- get_index_pars(chain, 1)
omg <- gaussian_process_model(tmp_pars, times)

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
