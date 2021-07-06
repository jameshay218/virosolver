#' Plot probability of infection
#' 
#' Plot probabilities of infection from compartmental model. Returns the probabilities and the plot.
#' 
#' @param chain A dataframe containing the MCMC samples 
#' @param nsamps Number of samples
#' @param INCIDENCE_FUNC A pointer to the Gaussian process model
#' @param solve_times Vector indicating the time over which the model is solved
#' @param obs_dat A dataframe containing observed Ct values and time of sample 
#' collection. NULL by default.
#' @param true_prob_infection A dataframe from simulated data with two columns, 
#' one for time and the other is the true probability of infection. 
#' NULL by default.
#' @param tshift Shift the solve times? Numeric, set to 0 by default
#' @param smooth Smooth the model estimates for plotting? FALSE by default.
#' 
#' @return Return a list containing three things: 1. A dataframe of model predictions containing
#' time, probability of infection, and sample number; 
#' 2. A dataframe containing the maximum posterior probability of infection and time;
#' 3. A ggplot showing the probabilities of infection
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family plots 
#' 
#' @examples
#'data(example_seir_incidence)
#'predictions <- plot_prob_infection(chain_comb, 
#'                                  nsamps=100, 
#'                                  INCIDENCE_FUNC=incidence_function,
#'                                  solve_times=0:max(ct_data_use$t),
#'                                  obs_dat=ct_data_use,
#'                                  true_prob_infection=example_seir_incidence)
#'p_incidence_prediction <- predictions$plot + scale_x_continuous(limits=c(0,200))
#'p_incidence_prediction
#' 
#' @export

plot_prob_infection <- function(chain,
                                nsamps,
                                INCIDENCE_FUNC,
                                solve_times,
                                obs_dat=NULL,
                                true_prob_infection=NULL,
                                tshift=0,
                                smooth=FALSE){
 
  ## Take n samples from the MCMC chain
   samps <- sample(unique(chain$sampno),nsamps)
  all_res <- NULL
  for(i in seq_along(samps)){
    samp <- samps[i]
    ## Return parameters for each sample according to sample number
    tmp_pars <- lazymcmc::get_index_pars(chain, samp)
    
  ## Solve compartmental model (SEIR or SEEIRR)
    prob_infection_tmp <- INCIDENCE_FUNC(tmp_pars, solve_times)
    
  ## Smooth out curve for plotting  
    if(smooth){
      prob_infection_tmp <- pmax(smooth.spline(prob_infection_tmp,spar=0.3)$y,0.0000001)
    }
    all_res[[i]] <- tibble(t=solve_times+tshift,prob_infection=prob_infection_tmp,sampno=i)
  }
  
  ## Combine results for all samples
  posterior_dat <- do.call("bind_rows",all_res)
  
  ## Return parameters with highest likelihood
  best_pars <- lazymcmc::get_best_pars(chain)
  
  ## Solve compartmental model using parameters with highest likelihood
  best_prob_infection <- INCIDENCE_FUNC(best_pars, solve_times)
  
  ## Smooth out curve for plotting  
  if(smooth){
    best_prob_infection <- pmax(smooth.spline(best_prob_infection)$y,0.0000001)
  }
  best_prob_dat <- tibble(t=solve_times+tshift,prob_infection=best_prob_infection,sampno="MAP")

  p1 <- ggplot(posterior_dat) +
    ## Plot probability of infection for each posterior sample
    geom_line(aes(x=t,y=prob_infection,group=sampno,col="Posterior draw"),size=0.1) +
    ## Plot one line for the MAP (maximum posterior probability of infection). This is the trajectory with the 
    ## highest likelihood.
    geom_line(data=best_prob_dat,aes(x=t,y=prob_infection,col="MAP"),size=0.5) +
    scale_y_continuous(expand=c(0,0))+
    xlab("Time") +
    ylab("Probability of infection") +
    theme_classic()

  ## Add line to the plot with the true probability of infection
  if(!is.null(true_prob_infection)){
    p1 <- p1 +
      geom_line(data=true_prob_infection,aes(x=t,y=prob_infection,col="Ground truth"),
                linetype="dashed",size=0.5)
  }

  ## Add vertical line with the sample date to the plot
  if(!is.null(obs_dat)){
    p1 <- p1 +
      geom_vline(data=data.frame(x=unique(obs_dat$t)),
                 aes(xintercept=x,linetype="Sample date"),
                 col="red",size=0.25)
  }
  
  ## Adding colors to the plot and modifying the legend
  ## Also returning the plot and a list containing probabilities of infection
  p1 <- p1 +
    scale_color_manual(values=c("Posterior draw"="gray50","MAP"="green","Ground truth"="blue")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position="bottom")
  return(list(predictions=posterior_dat, map_prediction=best_prob_dat, plot=p1))
}

#' Predicted distribution fits
#' 
#' Obtain predicted Ct distribution fits from model (posterior_dat)
#' 
#' @param chain A dataframe containing the MCMC samples 
#' @param MODEL_FUNC Function that expects a vector of model parameters with names 
#' corresponding to the parameter control table and returns a single log posterior probability
#' @param nsamps Number of samples. Defaults to 100.
#' 
#' @return Returns a dataframe containing the predictions from the posterior distribution.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family plots 
#' 
#' @examples
#' data(example_ct_data)
#' data(example_seir_partab)
#' 
#' \dontrun{
#'
#' MODEL_FUNC <- create_posterior_func(parTab=example_seir_partab,
#'                                      data=example_ct_data,
#'                                      PRIOR_FUNC=prior_func_seir,
#'                                      INCIDENCE_FUNC=incidence_function,
#'                                      use_pos=FALSE) 
#'                                      
#' posterior_dat <- predicted_distribution_fits(chain, MODEL_FUNC, nsamps=100)                                   
#' head(posterior_dat)
#' }
#' 
#' @export

predicted_distribution_fits <- function(chain, MODEL_FUNC, nsamps=100){
  
  ## Return model parameters with highest likelihood
  best_pars <- lazymcmc::get_best_pars(chain)

  ## Generate posterior draws for Ct distribution prediction
  ## Take n samples from the MCMC chain
  samps <- sample(unique(chain$sampno),nsamps)
  all_res <- NULL
  for(i in seq_along(samps)){
    samp <- samps[i]
    ## Return parameters for each sample according to sample number
    tmp_pars <- lazymcmc::get_index_pars(chain, samp)
    ## MODEL_FUNC output of create_posterior_func
    all_res[[i]] <- MODEL_FUNC(tmp_pars) %>% mutate(sampno=i)
  }
  
  ## Combine results for all samples
  posterior_dat <- do.call("bind_rows",all_res)
  return(posterior_dat)
}

#' Plot distribution fits 
#' 
#' Plot predicted Ct distribution fits from model
#' 
#' @param chain A dataframe containing the MCMC samples 
#' @param obs_dat A dataframe containing observed Ct values and time of sample 
#' collection. NULL by default.
#' @param MODEL_FUNC Function that expects a vector of model parameters with names 
#' corresponding to the parameter control table and returns a single log posterior probability
#' @param nsamps Number of samples. Defaults to 100.  
#' @param pos_only pos_only flag uses only Ct values below the limit of detection. Defaults to TRUE.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family plots 
#' 
#' @return Returns two stacked ggplots. 
#' 
#' @examples
#' \dontrun {
#' model_func_gp <- create_posterior_func(par_tab,example_ct_data,NULL,incidence_function,"model")
#' p_distribution_fit_gp <- plot_distribution_fits(chain_comb, example_ct_data, model_func_gp,100,pos_only=FALSE)
#' }
#' 
#' @export

plot_distribution_fits <- function(chain, obs_dat, MODEL_FUNC, nsamps=100, pos_only=TRUE){
  ## Pull out MAP parameters from MCMC chain
  best_pars <- get_best_pars(chain)
  
  ## Obtain predicted Ct distribution fits from model
  posterior_dat <- predicted_distribution_fits(chain, MODEL_FUNC, nsamps)
  
  ## Make sure only plotting detectable Ct values, and get label
  obs_dat1 <- obs_dat %>%
    filter(ct < best_pars["intercept"]) %>%
    mutate(obs_t=paste0("Sample day: ", t))

  ## Get number of observations per time point
  obs_tally <- obs_dat1 %>% group_by(t) %>% tally()

  ## Re-scale posterior densities to only detectable Ct distribution posterior densities
  total_density <- posterior_dat %>%
    filter(ct < best_pars["intercept"]) %>%
    group_by(t,sampno) %>%
    summarize(total_dens=sum(density)) %>% 
    left_join(obs_tally)

  ## Get expected number of observations for each Ct
  summary_posterior_dat <- posterior_dat %>%
    filter(ct < best_pars["intercept"]) %>%
    left_join(total_density) %>%
    filter(!is.na(n)) %>%
    group_by(t, sampno) %>%
    mutate(density=density/total_dens) %>%
    ungroup() %>%
    mutate(expectation=density*n) %>%
    ungroup() %>%
    ## Simulate observations
    mutate(sim_obs=rbinom(n(),n,density)) %>%
    group_by(ct, t)

  ## Summarise expected Ct values
  summary_expectation <- summary_posterior_dat %>%
    group_by(ct, t) %>%
    ## Obtain quantiles for expectations
    summarize(lower_expec=quantile(expectation,0.025),
              median_expec=quantile(expectation,0.5),
              upper_expec=quantile(expectation,0.975))

  ## Summarise observed Ct values
  summary_obs <- summary_posterior_dat %>%
    group_by(ct, t) %>%
    ## Obtain quantiles for observations
    summarize(lower_obs=quantile(sim_obs,0.025),
              median_obs=quantile(sim_obs,0.5),
              upper_obs=quantile(sim_obs,0.975))

  ## Plot model fit to detectable Ct distribution
  p1 <- ggplot(obs_dat1) +
    ## Plot histogram of observed Ct values
    geom_histogram(aes(x=ct,y=..count..),binwidth=1,fill="grey70",col="grey20",boundary=0) +
    ## 95% credible interval for observed Ct values
    geom_ribbon(data=summary_obs,aes(x=ct+0.5,ymin=lower_obs,ymax=upper_obs),fill="blue",alpha=0.25)+
    ## 95% credible interval for expected Ct values
    geom_ribbon(data=summary_expectation,aes(x=ct+0.5,ymin=lower_expec,ymax=upper_expec),fill="blue",alpha=0.5)+
    # median for expected Ct values
    geom_line(data=summary_expectation,aes(x=ct+0.5,y=median_expec),col="blue") +
    ## Reverse x-axis so Ct values decrease from high to low
    scale_x_continuous(trans="reverse",expand=c(0,0),limits=c(41,5),breaks=seq(0,40,by=5)) +
    coord_cartesian(xlim=c(0,39)) +
    coord_flip() +
    xlab("Ct value") +
    ylab("Count") +
    facet_wrap(~t,nrow=1,scales="free_x") +
    ## Style/Formatting 
    theme_classic() +
    theme(panel.grid.major = element_line(color="grey80",size=0.25),
          panel.grid.minor = element_line(color="grey80",size=0.25),
          panel.spacing = unit(0.5, "lines")) +
    ggtitle("Fit to detectable Ct distribution")

  ## Get predictions for probability undetectable
  summary_prop_detectable <- posterior_dat %>% filter(ct == best_pars["intercept"]) %>%
    group_by(t) %>%
    ## Density of prob undetectable is density of 1-prob of detectable
    mutate(density=1-density) %>%
    summarize(lower=quantile(density,0.025),
              median=quantile(density,0.5),
              upper=quantile(density,0.975))
  
  ## Prepare data to plot the proportion of detectable Ct values 
  p2 <-  ggplot(obs_dat %>%
                  group_by(t) %>%
                  mutate(is_detectable=ct < best_pars["intercept"]) %>%
                  summarize(prop_detectable=sum(is_detectable)/n()))
  
  ## pos_only flag uses only Ct values below the limit of detection
  if(!pos_only){
  p2 <- p2 +
    ## Plot proportion of detectable Ct values at t = 0.25 
    geom_point(aes(y=prop_detectable,x=0.25,col="Data"),size=3,shape=18)
  }
  
  ## Plot the proportion of detectable Ct values for positive and negative test results
  p2 <- p2 +
    geom_point(data=summary_prop_detectable,aes(x=0.75,y=median,col="Posterior median & 95% CI"),size=1) +
    ## Plot 95% CI 
    geom_errorbar(data=summary_prop_detectable,aes(x=0.75,ymin=lower,ymax=upper),
                  width=0.1, col="blue") +
    scale_y_continuous(limits=c(0,1)) +
    scale_x_continuous(limits=c(0,1)) +
    scale_color_manual(values=c("Data"="grey40",
                                "Posterior median & 95% CI"="blue",
                                "MAP"="green")) +
    guides(color=guide_legend(title=NULL)) +
    facet_wrap(~t,nrow=1) +
    ylab("Proportion detectable") +
    xlab("Sample time") +
    ## Style/Formatting
    theme_classic() +
    theme(panel.grid.major = element_line(color="grey80",size=0.25),
          panel.grid.minor = element_line(color="grey80",size=0.25),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(0.5, "lines"),
          axis.text.x=element_blank(),
          legend.position="none",
          axis.ticks.x=element_blank()) +
    ggtitle("Fit to proportion detectable")
    ## Lay out the plots using the R package "patchwork"
    return(p1/p2)
}

#' Plot posterior density
#'
#' Plot posterior density from model output 
#' 
#' @param chain 
#' @param var_name Character string denoting variable name from parameter table. 
#' @param parTab Dataframe specifying prior distribution for model parameter.
#' @param prior_mean Mean of the prior distribution for a particular parameter. 
#' @param prior_sd Standard deviation of the prior distribution for a particular parameter. 
#' @param real_data Plot true parameter values to simulate the data. Defaults to FALSE. 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family plot
#' 
#' @return ggplot
#' 
#' @examples
#' data(example_seir_partab)
#' chains <- load_mcmc_chains(location="mcmc_chains/readme_multiple_cross_section/",
#' parTab=par_tab,
#' burnin=mcmc_pars["adaptive_period"],
#' chainNo=FALSE,
#' unfixed=FALSE,
#' multi=FALSE)
#' 
#' var_name <- "tshift"
#' prior_mean <- subset(example_seir_partab, names == "tshift")$values 
#' prior_sd <- 0.2 
#' plot_posterior_density <- function(chain, var_name, parTab, prior_mean, prior_sd, real_data=FALSE)
#' 
#' @export

plot_posterior_density <- function(chain, var_name, parTab, prior_mean, prior_sd, real_data=FALSE){
  p <- ggplot(chain) +
    ## Plot density from MCMC chain (posterior distribution)
    geom_density(aes_string(var_name),fill="red",alpha=0.25) +
    ## Plot a normal function on top using the prior distribution 
    stat_function(data=data.frame(x=c(parTab[parTab$names == var_name,"lower_bound"],
                                      parTab[parTab$names == var_name,"upper_bound"])),
                  aes(x),col="blue",
                  fun = dnorm, n = 101, args = list(mean = prior_mean, sd = prior_sd)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Density") +
    theme_classic()
  ## Plot true parameter value used to simulate the data 
  if(!real_data){
   p <- p +
     geom_vline(xintercept=parTab[parTab$names == var_name,"values"], linetype="dashed")
  }
  return(p)
}

