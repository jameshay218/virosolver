#' Creates the posterior function to be used in the MCMC framework.
#' 
#' Creates a new function which calculates the posterior probability of a set 
#' of parameter values conditional on the Ct data.
#' 
#' @param parTab A vector of model parameters with names corresponding to the parameter 
#' control table.
#' @param data A tibble containing two columns: t for sample time and ct for the Ct values.
#' @param PRIOR_FUNC A function for the prior distribution. NULL by default.
#' @param INCIDENCE_FUNC A function that returns a vector of daily infection probabilities/
#' incidence. NULL by default.
#' @param solver_ver Can be set to 'model' or 'likelihood.' 'model' returns the predicted 
#' Ct distribution rather than the posterior probability. Set to 'likelihood' by default.
#' @param solve_likelihood Solves the likelihood. TRUE by default.
#' @param use_pos To fit the model using ALL PCR results (i.e. including samples testing
#'  negative), set this to FALSE. If the data only include positive Ct values (i.e. only 
#'  Ct values < the LOD), then this should be set to TRUE. FALSE by default.
#' 
#' @return Returns a single log posterior probability or the predicted Ct distribution.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family create posterior functions
#' 
#' @examples
#' incidence_function <- solveSEIRModel_lsoda_wrapper
#'data(example_seir_partab)
#'posterior_func <- create_posterior_func(parTab=example_seir_partab,
#'                                        data=example_ct_data,
#'                                        PRIOR_FUNC=prior_func_seir,
#'                                        INCIDENCE_FUNC=incidence_function,
#'                                        use_pos=FALSE) ## Important argument, see text
#'posterior_func(example_seir_partab$values)
#' 
#' @export

create_posterior_func <- function(parTab,
                                  data,
                                  PRIOR_FUNC=NULL,
                                  INCIDENCE_FUNC=NULL,
                                  solve_ver="likelihood",
                                  solve_likelihood=TRUE,
                                  use_pos=FALSE,
                                  ...) {
  par_names <- parTab$names
  pars <- parTab$values
  names(pars) <- par_names
  times <- 0:max(data$t)
  ages <- 1:max(data$t)
  obs_times <- unique(data$t)

  ## Pull out undetectable Ct values and count how many per observation time
  ## We only need get the likelihood of an undetectable Ct value once per time point,
  ## and then just have this contribute N times where N is the number of undetectable Cts
  ## at that time point
  data_use <- data
  undetectable_counts <- NULL
  if("intercept" %in% par_names){
    undetectable_tally <- data_use %>%
      filter(ct >= pars["intercept"]) %>%
      group_by(t) %>%
      tally()
    no_undetectable_times <- setdiff(obs_times, unique(undetectable_tally$t))
    no_undetectable_tally <- tibble(t=no_undetectable_times,n=0)
    undetectable_tally <- bind_rows(undetectable_tally, no_undetectable_tally) %>% arrange(t)
    undetectable_counts <- undetectable_tally$n
    data_use <- data_use %>% filter(ct < pars["intercept"])
  }


  ## Pull out data into a list for quicker indexing later on
  data_list <- NULL
  for(i in seq_along(obs_times)){
    data_list[[i]] <- data_use %>% filter(t == obs_times[i]) %>% pull(ct)
  }


  ## Function that returns either the likelihood or model predictions
  f <- function(pars){
    names(pars) <- par_names
    
    ## Get probability of infection
    prob_infection_tmp <- INCIDENCE_FUNC(pars, times)
    
    ## Returns the likelihood
    if(solve_ver == "likelihood"){
      lik <- 0
      if(solve_likelihood){
        lik <- sum(likelihood_cpp_wrapper(data_list, ages, obs_times,pars, prob_infection_tmp,
                                          use_pos,undetectable_counts))
      }

      ## Returns the sum of the likelihood and the prior distribution
      if(!is.null(PRIOR_FUNC)){
        prior <- PRIOR_FUNC(pars, ...)
        lik <- lik+prior
      }
      return(lik)
      
      ## Otherwise, return the model predictions
    } else {
      preds <- pred_dist_wrapper(seq(0,40,by=1),obs_times,ages,pars,prob_infection_tmp)
      return(preds)
    }
  }
  f
}

#' Creates the posterior function used in the MCMC framework for the SEEIRR model
#'
#' Creates a new function which calculates the posterior probability of a set 
#' of parameter values conditional on the Ct data using an SEEIRR model.
#' 
#' @param parTab A vector of model parameters with names corresponding to the parameter 
#' control table.
#' @param data A tibble containing three columns: date for the sample collection date, 
#' POS for the number of individuals who test positive, and NEG for the number of 
#' individuals who test negative. 
#' @param ts The times over which the SEEIRR model is solved. Used to subset the modeled
#' prevalence by the sample collection dates.
#' @param INCIDENCE_FUNC A function that returns a vector of daily infection probabilities/
#' incidence. Set to 'detectable_SEEIRRModel' by default.
#' @param PRIOR_FUNC A function for the prior distribution. NULL by default.
#' @param ver Can be set to 'model' or 'likelihood.' 'model' returns the prevalence
#' from the SEEIRR model rather than the posterior probability. Set to 'likelihood' by default.
#' 
#' @return Returns a single log posterior probability or the prevalence from the 
#' SEEIRR model.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family create posterior functions
#' 
#' @examples FIX ME
#'
#' @export

create_post_func_seeirr <- function(parTab, data, ts, INCIDENCE_FUNC=detectable_SEEIRRModel, PRIOR_FUNC=NULL,ver="likelihood"){
  par_names <- parTab$names
  test_times <- unique(data$date)
  observed_prev <- data$POS
  N_obs <- data$POS + data$NEG

  f <- function(pars){
    names(pars) <- par_names
    
    ## Solves the SEEIRR model for infection prevalence
    prev <- INCIDENCE_FUNC(pars, ts)
    
    if(ver == "likelihood"){
      
      ## Gets the prevalence for the times that correspond to the sample collection dates
      prev <- prev[which(ts %in% test_times)]
      
      ## Calculates the likelihood
      lik <- sum(dbinom(observed_prev, N_obs, prev,log=TRUE))
      
      ## Sums the likelihood and the prior
      if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
      
      return(lik)
      
      ## Otherwise, return the prevalence
    } else {
      return(prev)
    }
  }
  f
}

#' Creates the posterior function used in the MCMC framework for the combined SEEIRR model 
#' 
#' Creates a new function which calculates the posterior probability of a set 
#' of parameter values conditional on the Ct data using the combined SEEIRR model. This
#' model uses multiple cross-sectional samples.
#' 
#' #' @param parTab A vector of model parameters with names corresponding to the parameter 
#' control table.
#' @param data A tibble containing four columns: date for the sample collection date, 
#' POS for the number of individuals who test positive, NEG for the number of 
#' individuals who test negative, and location for location of sample collection.
#' @param ts The times over which the SEEIRR model is solved. Used to subset the modeled
#' prevalence by the sample collection dates.
#' @param INCIDENCE_FUNC A function that returns a vector of daily infection probabilities/
#' incidence. Set to 'detectable_SEEIRRModel' by default.
#' @param PRIOR_FUNC A function for the prior distribution. NULL by default.
#' @param ver Can be set to 'model' or 'likelihood.' 'model' returns the prevalence
#' from the SEEIRR model rather than the posterior probability. Set to 'likelihood' by default.
#' 
#' @return Returns a single log posterior probability or the prevalence from the combined
#' SEEIRR model.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family create posterior functions
#' 
#' @examples FIX ME
#' 
#' @export

create_post_func_seeirr_combined <- function(parTab, data, ts, INCIDENCE_FUNC=detectable_SEEIRRModel,
                                             PRIOR_FUNC=NULL,ver="likelihood"){
  par_names <- parTab$names
  unique_locs <- unique(data$location)
  f <- function(pars){
    names(pars) <- par_names
    all_lik <- numeric(length(unique_locs))
    all_prev <- NULL

    ## for each location, prepare a temporary dataset
    for(i in seq_along(unique_locs)){
      loc <- unique_locs[i]
      tmp_dat <- data %>% filter(location == loc)
      test_times <- unique(tmp_dat$date)
      observed_prev <- tmp_dat$POS
      N_obs <- tmp_dat$POS + tmp_dat$NEG

      ## Establish some temporary parameters
      tmp_pars <- c(pars[paste0("R0_",i)], pars[paste0("t0_",i)],
                    pars["infectious"],pars["latent"],pars["incubation"],pars["recovery"],pars["I0"])
      names(tmp_pars) <- c("R0","t0","infectious","latent","incubation","recovery","I0")

      ## Solve the SEEIRR model for infection prevalence
      prev <- INCIDENCE_FUNC(tmp_pars, ts)[1:length(ts)]
      
      ##  Get the prevalence for the times that correspond to the sample collection dates
      prev_subset <- prev[which(ts %in% test_times)]
      
      ## Calculate the likelihood
      lik <- sum(dbinom(observed_prev, N_obs, prev_subset,log=TRUE))
      
      ## Sum the likelihood and the prior
      if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
      
      ## Put the likelihood into a vector
      all_lik[i] <- lik
      
      ## Store the prevalence
      if(ver == "model"){
        all_prev[[i]] <- tibble(prev=prev,loc=loc,t=ts)
      }
    }
    
    ## Returns the total likelihood across all locations
    if(ver == "likelihood"){
      return(sum(all_lik))
    }
    
    ## Combines the prevalence across all locations
    if(ver == "model") {
      all_prev <- do.call("bind_rows", all_prev)
      return(all_prev)
    }
    
    return(all_lik)
  }
  f
}

#' Creates the posterior function used in the MCMC framework for detectable Ct values
#'
#' Creates a new function which calculates the posterior probability of a set 
#' of parameter values conditional on data that include only individuals with 
#' detectable PCR values on the sample collection day.
#' 
#' @param parTab A vector of model parameters with names corresponding to the parameter 
#' control table.
#' @param data A tibble containing two columns: t for sample time and ct for the Ct values.
#' @param PRIOR_FUNC A function for the prior distribution. NULL by default.
#' @param INCIDENCE_FUNC A function that returns a vector of daily infection probabilities/
#' incidence. NULL by default.
#' @param solve_likelihood Solves the likelihood. TRUE by default.
#'
#'@return Returns a single log posterior probability.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family create posterior functions
#' 
#' @examples FIX ME
#'
#' @export

create_posterior_func_detectable <- function(parTab,
                                  data,
                                  PRIOR_FUNC=NULL,
                                  INCIDENCE_FUNC=NULL,
                                  solve_likelihood=TRUE,
                                  ...) {
  par_names <- parTab$names
  pars <- parTab$values
  names(pars) <- par_names
  times <- 0:max(data$t)
  ages <- 1:max(data$t)
  obs_times <- unique(data$t)

  f <- function(pars){
    names(pars) <- par_names
    
    ## Get the probability of infection
    prob_infection_tmp <- INCIDENCE_FUNC(pars, times)
    
    lik <- 0
    
    ## Returns the likelihood
    if(solve_likelihood){
      lik <- sum(likelihood_detectable(data, ages, pars, prob_infection_tmp))
    }

    ## Returns the sum of the likelihood and the prior distribution
    if(!is.null(PRIOR_FUNC)){
      prior <- PRIOR_FUNC(pars, ...)
      lik <- lik+prior
    }
    
    return(lik)
  }
  f
}
