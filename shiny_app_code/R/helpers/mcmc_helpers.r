create_prior_func_seir <- function(seir_partab) {
  means <- seir_partab$values
  names(means) <- seir_partab$names
  ## Set standard deviations of prior distribution
  sds_seir <- c("obs_sd"=0.5,"viral_peak"=2,
                "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
                "prob_detect"=0.03,
                "incubation"=0.25, "infectious"=0.5)
  
  ## Define a function that returns the log prior probability for a given vector of parameter
  ## values in `pars`, given the prior means and standard deviations set above.
  prior_func_seir <- function(pars,...){
    ## Ct model priors
    obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_seir["obs_sd"],log=TRUE)
    viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_seir["viral_peak"],log=TRUE)
    wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_seir["wane_rate2"],log=TRUE)
    tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_seir["t_switch"],log=TRUE)
    level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_seir["level_switch"],log=TRUE)
    ## Beta prior on the prob_detect parameter to ensure between 0 and 1
    beta1_mean <- means[which(names(means) == "prob_detect")]
    beta1_sd <- sds_seir["prob_detect"]
    beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
    beta_beta <- beta_alpha*(1/beta1_mean - 1)
    beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
    
    ## SEIR model priors
    incu_prior <- dlnorm(pars["incubation"],log(means[which(names(means) == "incubation")]), sds_seir["incubation"], TRUE)
    infectious_prior <- dlnorm(pars["infectious"],log(means[which(names(means) == "infectious")]),sds_seir["infectious"],TRUE)
    
    ## Sum up
    obs_sd_prior + viral_peak_prior + 
      wane_2_prior + tswitch_prior + level_prior + beta_prior +
      incu_prior + infectious_prior
  }
}

single_xs <- function(time, ct_threshold=40, ct_data, 
                      seir_partab,
                      singlemodal=TRUE,
                      iterations=200000, #Actual default:200000
                      adaptive_period=100000, #Actual default:100000
                      nchains=3) {
  print("Constructing prior SEIR function...")
  prior_func_seir <- create_prior_func_seir(seir_partab)
  print("Setting incidence function and creating posterior function...")
  incidence_function <- solveSEIRModel_lsoda_wrapper
  ## Create the posterior function used in the MCMC framework
  posterior_func <- create_posterior_func(parTab=seir_partab,
                                          data=ct_data,
                                          PRIOR_FUNC=prior_func_seir,
                                          INCIDENCE_FUNC=incidence_function,
                                          use_pos=FALSE)
  print("Filtering ct data...")
  ct_data_use <- ct_data %>% filter(t==time)
  p_ct_use <- ggplot(ct_data_use %>% filter(ct < ct_threshold)) + 
    geom_histogram(aes(x=ct)) + theme_bw()
  vals <- list(p_ct_use)
  print("Generating starting parameter table...")
  start_tab <- virosolver::generate_viable_start_pars(seir_partab,ct_data_use,
                                          create_posterior_func,
                                          incidence_function,
                                          prior_func_seir)
  covMat <- diag(nrow(start_tab))
  mvrPars <- list(covMat,2.38/sqrt(nrow(start_tab[start_tab$fixed==0,])),w=0.8)
  mcmc_pars <- c("iterations"=iterations,"popt"=0.234,"opt_freq"=2000,
                 "thin"=1000,"adaptive_period"=adaptive_period,"save_block"=100)
  dir.create("mcmc_chains/readme_single_cross_section",recursive=TRUE)
  
  ##################################
  ## RUN THE MCMC FRAMEWORK'
  browser()
  #stdout <- vector('character')
  #con <- textConnection('stdout', 'wr', local = TRUE)
  #sink(con)
  
  res <- foreach(chain_no=1:nchains,.packages = c("virosolver","lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
    run_MCMC(parTab=start_tab,
                      data=ct_data_use,
                      INCIDENCE_FUNC=incidence_function,
                      PRIOR_FUNC=prior_func_seir,
                      mcmcPars=mcmc_pars,
                      filename=paste0("mcmc_chains/readme_single_cross_section/readme_seir_",chain_no),
                      CREATE_POSTERIOR_FUNC=create_posterior_func,
                      mvrPars=mvrPars,
                      use_pos=FALSE) ## Important argument
    #sink()
    #close(con)
    chains <- load_mcmc_chains(location="mcmc_chains/readme_single_cross_section",
                               parTab=start_tab,
                               burnin=mcmc_pars["adaptive_period"],
                               chainNo=TRUE,
                               unfixed=TRUE,
                               multi=TRUE)
    
    chains_melted <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))
    ## Look at trace plots
    print("Generating trace plots...")
    p_trace <- ggplot(chains_melted) + 
      geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
      facet_wrap(~name,scales="free_y") + 
      scale_color_viridis_d(name="Chain") + 
      theme_bw() +
      xlab("Iteration") +
      ylab("Value")
    vals <- append(vals, p_trace)
    ## Load in MCMC chains again, but this time read in the fixed parameters too 
    ## to ensure that the posterior draws are compatible with the model functions
    chains <- load_mcmc_chains(location="mcmc_chains/readme_single_cross_section",
                               parTab=start_tab,
                               burnin=mcmc_pars["adaptive_period"],
                               chainNo=FALSE,
                               unfixed=FALSE,
                               multi=TRUE)
    ## Do some reshaping to allow correct subsampling (we need each sampno to correspond to one unique posterior draw)
    chain_comb <- chains$chain %>% as_tibble() %>% mutate(sampno=1:n()) %>% as.data.frame()
    
    ## Load in true incidence curve to compare to our prediction
    data(example_seir_incidence)
    print("Generating predictions...")
    predictions <- plot_prob_infection(chain_comb, 
                                       nsamps=100, 
                                       INCIDENCE_FUNC=incidence_function,
                                       solve_times=0:max(ct_data_use$t),
                                       obs_dat=ct_data_use,
                                       true_prob_infection=example_seir_incidence)
    p_incidence_prediction <- predictions$plot + scale_x_continuous(limits=c(0,200))
    vals <- append(vals, p_incidence_prediction)
    print("Plotting model-predicted Ct distribution against sample data...")
    ## Use create_posterior_func to return the predicted Ct distribution rather than the posterior probability
    model_func <- create_posterior_func(example_seir_partab,ct_data_use,NULL,incidence_function,"model")
    ## Pass model_func to a plotting function to observe predicted Ct distribution against data
    p_distribution_fit <- plot_distribution_fits(chain_comb, ct_data_use, model_func,100,pos_only=FALSE)
    ## Joining, by = "t"
    ## Joining, by = c("t", "sampno")
    ## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
    vals <- append(vals, p_distribution_fit)
    }
  
  
}