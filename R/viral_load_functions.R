#' This function returns modal Ct values or viral load values, if 
#' convert_vl flag is set to TRUE. 
#' 
#' @param pars Vector of model parameter values. 
#' @param obs_t Vector of times since infection. Referred to as "ages."
#' @param convert_vl Boolean flag used to indicate if Ct values should
#' be converted to viral load values. Set to FALSE by default. 
#' @param infection_time  Time of infection. Set to 0 by default. 
#' 
#' @return Returns modal Ct or viral load values. 
#' 
#' @author James Hay, \email{jameshay218@@gmail.com}
#' @family viral load functions
#' 
#' @export
viral_load_func <- function(pars, obs_t, convert_vl=FALSE, infection_time=0){
  ## Days post infection until growth
  tshift <- pars["tshift"] + infection_time
  ## Days post growth to peak
  desired_mode <- pars["desired_mode"]
  ## Days post peak to switch point
  t_switch <- pars["t_switch"] 
 
  ## Peak Ct value
  height <- pars["viral_peak"]
  ## Ct value at the switch
  level_switch <- pars["level_switch"] 

  ## y-axis shift
  true_0 <- pars["true_0"] 
  ct_intercept <- pars["intercept"] 
  wane_rate <- (height-level_switch)/t_switch
  wane_rate2 <- (level_switch-true_0)/pars["wane_rate2"]
  growth_rate <- (height-true_0)/desired_mode

  ## Categorizes observed Cts into time periods
  t_period1 <- obs_t <= tshift
  t_period2 <- obs_t > tshift & obs_t <= (desired_mode + tshift)
  t_period3 <- obs_t > (desired_mode + tshift) & obs_t <= (desired_mode + tshift + t_switch)
  t_period4 <- obs_t > (desired_mode + tshift + t_switch)
  y <- numeric(length(obs_t))
  
  ## Calculates and assigns modal Ct based on time period 
  y[t_period1] <- true_0
  y[t_period2] <- growth_rate * (obs_t[t_period2] - tshift) + true_0
  y[t_period3] <- height - wane_rate * (obs_t[t_period3] - (desired_mode+tshift))
  y[t_period4] <- level_switch - wane_rate2*(obs_t[t_period4] - (desired_mode + tshift + t_switch))
  ct <- y
  ## Converts Ct values to viral load
  if(convert_vl){
    vl <- ((ct_intercept-ct)/log2(10)) + pars["LOD"]
    y <- vl
  }
  y
}

#' This function calculates predicted Ct value densities for each observed
#' time and returns a tibble of Ct values, predicted Ct value densities, 
#' and observed times. 
#' 
#' @param test_cts Vector of Ct values.
#' @param obs_times Vector of times at which samples were collected (cross sections).
#' @param ages Vector of ages (time since infection).
#' @param pars Model parameters.
#' @param prob_infection Vector of the probability of infection.
#' @param symptom_surveillance Boolean if TRUE, then takes arguments from pars and generates the Ct distribution assuming symptom-based surveillance
#' @param sampling_dist character, either "gamma" or "uniform" giving the distribution used for the sampling delay distribution
#' 
#' @return Tibble containing Ct values, densities, and times 
#' 
#' @author James Hay, \email{jameshay218@@gmail.com}
#' @family viral load functions
#' 
#' @export
pred_dist_wrapper <- function(test_cts, obs_times, ages, pars, prob_infection, symptom_surveillance=FALSE, sampling_dist="gamma"){
  max_age <- length(ages)
  
  if(symptom_surveillance){
    test_cts <- test_cts[test_cts < pars["intercept"]]
    max_age <- pars["max_incu_period"] + pars["max_sampling_delay"] + 1
  }
  
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], max_age) ## Up until t_switch, full standard deviation
  
  ## Prior to t_switch, sd=1
  unmod_vec <- 1:min(t_switch,max_age)
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, variance about modal Ct trajectory decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])

  if(sampling_dist == "gamma"){
    sampling_dist_int <- 1
  } else {
    sampling_dist_int <- 2
  }
  
  comb_dat <- NULL
  for(obs_time in obs_times){
    ## Restrict ages to times occurring before 
    ## time of sample collection (single cross section)
    ages1 <- ages[(obs_time - ages) > 0]
    ## Returns the full probability density distribution to simulate from
    if(!symptom_surveillance){
      densities <- pred_dist_cpp(test_cts, ages1, obs_time, pars, prob_infection,sd_mod) 
    } else {
      densities <- pred_dist_cpp_symptoms(test_cts, pars["max_incu_period"],pars["max_sampling_delay"], obs_time, pars, prob_infection,sd_mod,sampling_dist_int)
    }
    comb_dat[[obs_time]] <- tibble(ct=test_cts,density=densities, t=obs_time)
  }
  comb_dat <- do.call("bind_rows",comb_dat)
  comb_dat
}
