#' Wrapper function for the likelihood code, which is written in C. 
#' 
#' This function returns the likelihood of observing a particular set of 
#' Ct values conditional on the underlying model parameters.
#' 
#' @param obs_dat Tibble of observed data. Contains two columns, t for time of sample
#' collection and ct for the Ct values.
#' @param ages Vector of times since infection.
#' @param times Vector of sample collection times in the data.
#' @param pars Model parameters.
#' @param prob_infection Probability of infection.
#' @param pos_only only use detectable Ct values? FALSE by default.
#' @param undetectable_counts Number of undetectable Ct samples. NULL by default.
#' 
#' @return Returns the likelihood.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples FIX ME
#' 
#' @export

likelihood_cpp_wrapper <- function(obs_dat, ages, times,
                                   pars, prob_infection,
                                   pos_only=FALSE,
                                   undetectable_counts=NULL){
  
  ## If flag to only use detectable Ct values or not (use different likelihood function as needed)
  liks_tj <- 0
  if(pos_only) {
    use_func <- likelihood_pos_only_cpp
  } else {
    use_func <- likelihood_cpp
  }

  ## Because we have a different standard deviation for different times
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], max(ages))

  ## Prior to t_switch, full standard deviation
  ## Make sure we don't go past maximum of ages vector
  unmod_vec <- 1:min(t_switch,max(ages)+1)
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  ## The rest are at sd_mod
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
  

  ## For each sampling times, calculate likelihood for Ct values measured at that time
  for(i in seq_along(times)){
    ## Only solve back until the earliest possible infection time
    ages1 <- ages[(times[i] - ages) > 0]
    sd_mod1 <- sd_mod[(times[i] - ages) > 0]
    obs1 <- obs_dat[[i]]

    ## Undetectable Cts can be bucketed, as same likelihood for each
    undetectable_lik <- 0
    if(!is.null(undetectable_counts) & !pos_only) {
      undetectable_lik <-  use_func(pars["intercept"], times[i], ages1,
                                    pars, prob_infection,
                                    sd_mod)*undetectable_counts[i]
    }

    ## Sum the likelihoods across all times
    liks_tj <- liks_tj +
      sum(use_func(obs1, times[i], ages1, pars, prob_infection,sd_mod)) + ## Detectable Cts
      undetectable_lik
  }
  liks_tj
}

#' Probability curve of detectable Ct values.
#' 
#' This function takes the model parameters and returns a gamma-distributed curve of 
#' the probabilities of detectable viral loads with parameters pars and ages.
#' 
#' @param pars Model parameters.
#' @param ages Vector of times since infection.
#' 
#' @return Returns a gamma-distributed curve of the probabilities of detectable 
#' viral loads.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples FIX ME
#' 
#' @export

prob_detectable_curve <- function(pars, ages){
  
  ## Extract parameter values
  mean <- pars["mean"]
  var <- pars["var"]

  ## Convert parameters for gamma distribution
  scale <- var/mean
  shape <- mean/scale

  ## Call the gamma distribution using model parameters
  probs <- dgamma(ages,shape=shape,scale=scale)*pars["c"]
  
  ## Make sure probabilities do not exceed 1
  probs[probs > 1] <- 1

  return(probs)

}

#' Likelihood for using proportion detectable only
#'
#' This function returns the likelihood of observing a particular set of Ct 
#' values conditional on the underlying model parameters using only the proportion
#' of Ct samples that are detectable.
#' 
#' @param obs_dat Tibble of observed data. Contains two columns, t for time of sample
#' collection and ct for the Ct values.
#' @param ages Vector of times since infection.
#' @param pars Model parameters.
#' @param prob_infection Probability of infection.
#' 
#' @return Returns the likelihood.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples FIX ME
#'
#' @export

likelihood_detectable <- function(obs_dat, ages, pars, prob_infection){
  liks <- 0
  ## For each sampling times, calculate likelihood for Ct values measured at that time
  for(i in 1:nrow(obs_dat)){
    ## Only solve back until the earliest possible infection time
    obs_time <- obs_dat$t[i]
    ## 'ages' indicates a range of times since infection
    ages1 <- ages[(obs_time - ages) > 0]
    n_pos <- obs_dat$pos[i]
    n_neg <- obs_dat$neg[i]
    ## Sum of detectable infections
    sum_detect_inf <- sum(vapply(ages1,
                                 FUN=function(a) prob_infection[obs_time-a]*prob_detectable_curve(pars, a),
                                 FUN.VALUE=numeric(1)))
    ## Likelihood of detectable 
    liks_detect <- n_pos*log(sum_detect_inf)
    ## Likelihood of undetectable
    liks_undetect <- n_neg*log(1-sum_detect_inf)
    ## Sum likelihoods of detectable and undetectable
    liks <- liks +  liks_detect+liks_undetect

  }
  liks
}


#' Function that gives the probability of observing Ct value x given age 
#' (days after infection) a and the viral kinetics curve.
#' 
#' This is the Gumbel probability density function normalized to the observable values.
#' 
#' @param x Observed Ct value.
#' @param a Days after infection.
#' @param pars Model parameters.
#' @param viral_loads Observed viral loads.
#' @param sd_mod Standard deviation of the viral kinetics curve.
#' 
#' @return Returns a vector of probabilities. 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples FIX ME
#' 
#' 
#' @export

p_a <- function(x, a, pars, viral_loads, sd_mod) {
  viral_load_sd <- pars["obs_sd"]*sd_mod[a]
  ## The intercept is the maximum number of cycles run on the PCR machine
  LOD <- pars["intercept"]
  
  renormalize <- extraDistr::pgumbel(LOD, viral_loads[a],viral_load_sd,lower.tail=TRUE) -
    extraDistr::pgumbel(0, viral_loads[a],viral_load_sd,lower.tail=TRUE)
  
  probs <- extraDistr::dgumbel(x, mu=viral_loads[a], sigma=viral_load_sd, log=FALSE)/renormalize
  
  probs
}

#' Proportion still detectable (single time)
#' 
#' Gives the proportion of infections still detectable for a single time since infection (a)
#' 
#' @param a Days after infection.
#' @param pars Model parameters.
#' @param viral_loads Observed viral loads.
#' @param sd_mod Standard deviation of the viral kinetics curve.
#' 
#' @return Returns a vector of probabilities. 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples FIX ME
#' 
#' @export

prop_detectable_single <- function(a, pars, viral_loads, sd_mod){
  viral_load_sd <- pars["obs_sd"]*sd_mod[a]
  LOD <- pars["intercept"]
  additional_prob <- 1
  ## How many days spent in detectability loss phase?
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]

  days_potential_loss <- a - t_switch
  
  ## Probability of additional decline in detectability
  if(days_potential_loss >= 0){
    additional_prob <- (1-pars["prob_detect"]*pars["t_unit"])^days_potential_loss
  }
  
  main_probs <- extraDistr::pgumbel(LOD,mu=viral_loads[a],sigma=viral_load_sd, lower.tail=TRUE, log.p=FALSE)

  main_probs * additional_prob
}

#' Proportion still detectable (range of times)
#' 
#' Gives the proportion of infections still detectable over a range of times 
#' since infection (referred to as 'ages')
#' 
#' @param ages Vector of times since infection.
#' @param pars Model parameters.
#' @param viral_loads Observed viral loads.
#' 
#' @return Returns a vector of probabilities. 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples 
#' data(example_gp_partab)
#' pars <- example_gp_partab$values
#' names(pars) <- example_gp_partab$names
#' test_ages <- seq(1,50,by=1)
#' cts <- viral_load_func(pars, test_ages)
#' prop_detect <- prop_detectable(test_ages,pars, cts)
#' @export

prop_detectable <- function(ages, pars, viral_loads){
  
  ## Because we have a different standard deviation for different times
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], max(ages))

  ## Prior to t_switch, full standard deviation
  ## Make sure we don't go past maximum of ages vector
  unmod_vec <- 1:min(t_switch,max(ages)+1)
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  ## The rest are at sd_mod
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])

  ## Use the function for the proportion still detectable for a single time point
  vapply(ages, function(x){
    prop_detectable_single(x, pars, viral_loads, sd_mod)
  },FUN.VALUE=numeric(1))
}

#' Likelihood (R version)
#' 
#' This is the likelihood function written in R code. It returns the likelihood 
#' of observing a particular set of Ct values conditional on the underlying model 
#' parameters. Despite being slower than the corresponding likelihood function 
#' written in C, it more closely matches the written equations in the paper 
#' and is retained in the package for didactic purposes. 
#'
#' @param obs_dat Tibble of observed data. Contains two columns, t for time of sample
#' collection and ct for the Ct values.
#' @param ages Vector of times since infection.
#' @param pars Model parameters.
#' @param prob_infection Probability of infection.
#'
#' @return Returns the likelihood.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples FIX ME
#'
#' @export

likelihood_R <- function(obs_dat, ages, pars, prob_infection){
  LOD <- pars["intercept"]

  ## Because we have a different standard deviation for different times
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], length(ages)) ## Up until t_switch, full standard deviation

  ## Prior to t_switch, 1
  unmod_vec <- 1:min(t_switch,length(ages))
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  ## The rest are at sd_mod
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
 
  ## Get the modal Ct value
  viral_loads <- viral_load_func(pars, ages)
  
  prob_detectable_dat <- sapply(ages, function(a) prop_detectable_cpp(a, viral_loads[a],
                                                                      pars["obs_sd"]*sd_mod[a], pars["intercept"],
                                                                      t_switch, pars["prob_detect"]))

  times <- unique(obs_dat$t)
  lik_tj <- 0
  for(obs_time in times){
    ages1 <- ages[(obs_time - ages) > 0]
    obs1 <- obs_dat %>% filter(t == obs_time) %>% pull(ct)
    
    ## Probability undetectable
    undetectable_prob <- 1 - sum(sapply(ages1, function(a) prob_detectable_dat[a]*prob_infection[obs_time-a]))

    ## Probability detectable
    detectable_prob <- function(X_i){
      sum(vapply(ages1,
                 FUN=function(a) p_a(X_i, a, pars, viral_loads,sd_mod)*prob_infection[obs_time-a]*prob_detectable_dat[a],
                 FUN.VALUE=numeric(1)))
    }
    
    ## Likelihood
    lik_tj <- lik_tj + log(sapply(obs1,function(x_i) (x_i >= LOD)*undetectable_prob + (x_i < LOD)*detectable_prob(x_i)))
  }
  lik_tj
}


#' Likelihood for using proportion detectable only (R version)
#'
#' This function is written in R. It returns the likelihood of observing a particular set of Ct 
#' values conditional on the underlying model parameters using only the proportion
#' of Ct samples that are detectable.Despite being slower than the corresponding likelihood 
#' function written in C, it more closely matches the written equations in the paper and 
#' is retained for didactic purposes. 
#' 
#' @param obs_dat Tibble of observed data. Contains two columns, t for time of sample
#' collection and ct for the Ct values.
#' @param ages Vector of times since infection.
#' @param pars Model parameters.
#' @param prob_infection Probability of infection.
#'
#' @return Returns the likelihood.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family likelihood functions
#' 
#' @examples FIX ME
#' 
#' @export

likelihood_pos_R <- function(obs_dat, ages, pars, prob_infection){
  ## The intercept is the maximum number of cycles run on the PCR machine
  LOD <- pars["intercept"]
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]

  ## Get the modal Ct value
  viral_loads <- viral_load_func(pars, ages)
  
  prob_detectable_dat <- sapply(ages, function(a) prop_detectable_cpp(a, viral_loads[a],
                                                                      pars["obs_sd"], pars["intercept"],
                                                                      t_switch, pars["prob_detect"]))
  times <- unique(obs_dat$t)
  lik_tj <- 0
  for(obs_time in times){
    ages1 <- ages[(obs_time - ages) > 0]
    obs1 <- obs_dat %>% filter(t == obs_time) %>% pull(ct)
    
    ## Probability that an individual has a detectable Ct value and is infected
    prob_detectable_and_infected <- sum(sapply(ages1, function(a) prob_detectable_dat[a]*prob_infection[obs_time-a]))

    ## Probability of detectable Ct value
    detectable_prob <- function(X_i){
      sum(vapply(ages1,
                 FUN=function(a) p_a(X_i, a, pars, viral_loads,sd_mod)*prob_infection[obs_time-a]*prob_detectable_dat[a],
                 FUN.VALUE=numeric(1)))
    }
    ## Likelihood
    lik_tj <- lik_tj + log(sapply(obs1,function(x_i) detectable_prob(x_i)/prob_detectable_and_infected))
  }
  lik_tj
}


