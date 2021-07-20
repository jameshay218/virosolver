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
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family viral load functions
#' 
#' @examples FIX ME
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
  wane_rate2 <- (level_switch-ct_intercept)/pars["wane_rate2"]
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
#' 
#' @return Tibble containing Ct values, densities, and times 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family viral load functions
#' 
#' @examples FIX ME
#' 
#' @export
pred_dist_wrapper <- function(test_cts, obs_times, ages, pars, prob_infection, symptom_surveillance=FALSE){
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], length(ages)) ## Up until t_switch, full standard deviation

  ## Prior to t_switch, sd=1
  unmod_vec <- 1:min(t_switch,length(ages))
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, variance about modal Ct trajectory decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])

  comb_dat <- NULL
  for(obs_time in obs_times){
    ## Restrict ages to times occurring before 
    ## time of sample collection (single cross section)
    ages1 <- ages[(obs_time - ages) > 0]
<<<<<<< HEAD
    ## Returns the full probability density distribution to simulate from
    densities <- pred_dist_cpp(test_cts, ages1, obs_time, pars, prob_infection,sd_mod)  
=======
    if(!symptom_surveillance){
      densities <- pred_dist_cpp(test_cts, ages1, obs_time, pars, prob_infection,sd_mod) ## FIX ME: is pred_dist_cpp still relevant?  
    } else {
      densities <- pred_dist_cpp_symptoms(test_cts, pars["max_incu_period"],pars["max_sampling_delay"], obs_time, pars, prob_infection,sd_mod)
    }
>>>>>>> f5a70dd60ef03d6aaef9b88b0dfec542fe5b6d5f
    comb_dat[[obs_time]] <- tibble(ct=test_cts,density=densities, t=obs_time)
  }
  comb_dat <- do.call("bind_rows",comb_dat)
  comb_dat

}

#' This function simulates N Ct values per age in ages 
#' and returns a dataframe of Ct values and ages. 
#'  
#' @param ages Vector of ages (days since infection).
#' @param kinetics_pars Vector of named parameters for the viral kinetics model.
#' @param N Number of observations to randomly generate per age;
#' Set to 100 by default.
#' 
#' @return Dataframe of simulated Ct values and corresponding ages.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family viral load functions
#' 
#' @examples FIX ME
#' 
#' @export
simulate_viral_loads_example <- function(ages, kinetics_pars,N=100){
  t_switch <-  kinetics_pars["t_switch"] + kinetics_pars["desired_mode"] + kinetics_pars["tshift"]
  sd_mod <- rep(kinetics_pars["sd_mod"], max(ages))
  unmod_vec <- 1:min(t_switch,max(ages))
  sd_mod[unmod_vec] <- 1
  decrease_vec <- (t_switch+1):(t_switch+kinetics_pars["sd_mod_wane"])
  ## For the next sd_mod_wane days, variance about Ct trajectory decreases linearly
  sd_mod[decrease_vec] <- 1 - ((1-kinetics_pars["sd_mod"])/kinetics_pars["sd_mod_wane"])*seq_len(kinetics_pars["sd_mod_wane"])
  
  sim_dat <- matrix(ncol=N,nrow=length(ages))
  for(age in ages){
    ## For each value we're going to simulate, pre-determine if it will still be detectable or not
    detectable_statuses <- rnbinom(N, 1, prob=kinetics_pars["prob_detect"]) + 
      kinetics_pars["tshift"] + kinetics_pars["desired_mode"] + kinetics_pars["t_switch"]
    cts <- rep(virosolver::viral_load_func(kinetics_pars, age, FALSE, 0),N)
    cts[detectable_statuses <= age] <- 1000
    sd_used <- kinetics_pars["obs_sd"]*sd_mod[age]
    ## Generate N observations of Ct values from gumbel distribution for a specified mode
    ct_obs_sim <- extraDistr::rgumbel(N, cts, sd_used)
    ## Set Ct values greater than intercept to intercept value
    ct_obs_sim <- pmin(ct_obs_sim, kinetics_pars["intercept"])
    sim_dat[age,] <- ct_obs_sim
  }
  sim_dat <- data.frame(sim_dat)
  colnames(sim_dat) <- 1:N
  sim_dat$time_since_infection <- ages
  sim_dat <- sim_dat %>% pivot_longer(-time_since_infection) %>% rename(ct=value,age=time_since_infection,i=name)
  return(sim_dat)
}
