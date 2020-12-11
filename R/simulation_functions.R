#' @export
simulate_viral_loads <- function(infection_times,
                                 solve_times,
                                 kinetics_pars,
                                 vl_func_use=viral_load_func,
                                 convert_vl=FALSE,
                                 add_noise=NULL){
  viral_loads <- matrix(-100, nrow=length(infection_times), ncol=length(solve_times))
  n <- length(infection_times)


  ## Control for changing standard deviation
  test_ages <- 1:1000
  t_switch <-  kinetics_pars["t_switch"] + kinetics_pars["desired_mode"] + kinetics_pars["tshift"]
  sd_mod <- rep(kinetics_pars["sd_mod"], max(test_ages))
  unmod_vec <- 1:min(t_switch,max(test_ages))
  sd_mod[unmod_vec] <- 1
  decrease_vec <- (t_switch+1):(t_switch+kinetics_pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-kinetics_pars["sd_mod"])/kinetics_pars["sd_mod_wane"])*seq_len(kinetics_pars["sd_mod_wane"])


  ## Pre-compute negative binomial draws for loss in detectability
  if(additional_detect_process){
    ## How many days do you remain detectable after hitting the switch point?
    days_still_detectable <- rnbinom(n, size=1,prob=kinetics_pars["prob_detect"])
  }

  for(i in seq_along(infection_times)){
    if(i %% 1000 == 0) message(i)
    if(infection_times[i] > 0){
      pars <- kinetics_pars
      mod_probs <- rep(1, length(solve_times))
      ct <- vl_func_use(pars, times, FALSE, infection_times[i])
      ## Additional way that individuals can become undetectable
      if(additional_detect_process){
        ## How long to wait until undetectable?
        mod_probs[which(times >= kinetics_pars["tshift"] +
                          kinetics_pars["desired_mode"] +
                          kinetics_pars["t_switch"] +
                          infection_times[i] +
                          days_still_detectable[i])] <- 0
      }
      ct[ct > kinetics_pars["true_0"]] <- kinetics_pars["true_0"]
      ct[which(times < infection_times[i])] <- 1000
      ct <- ct * mod_probs
      ct[which(mod_probs == 0)] <- 1000
      viral_loads[i,] <- ct
    }
  }
  colnames(viral_loads) <- solve_times

  if(convert_vl) {
    viral_loads <- ((kinetics_pars["intercept"]-viral_loads)/log2(10)) + kinetics_pars["LOD"]
    true_viral_loads <- viral_loads
    true_viral_loads[true_viral_loads < kinetics_pars["LOD"]] <- kinetics_pars["LOD"]
  } else {
    true_viral_loads <- viral_loads
    true_viral_loads[true_viral_loads > kinetics_pars["intercept"]] <- kinetics_pars["intercept"]
  }
  if(!is.null(add_noise)){
    observed_viral_loads <- t(apply(viral_loads, 1, function(vl) add_noise(length(vl),vl, kinetics_pars["obs_sd"])))
    colnames(observed_viral_loads) <- solve_times
    if(convert_ct) {
      observed_viral_loads[observed_viral_loads < kinetics_pars["LOD"]] <- kinetics_pars["LOD"]
    } else {
      observed_viral_loads[observed_viral_loads > kinetics_pars["intercept"]] <- kinetics_pars["intercept"]
    }
  } else {
    observed_viral_loads <- true_viral_loads
  }
  true_viral_loads_melted <- reshape2::melt(true_viral_loads)
  colnames(true_viral_loads_melted) <- c("i","t","true")
  obs_viral_loads_melted <- reshape2::melt(observed_viral_loads)
  colnames(obs_viral_loads_melted) <- c("i","t","obs")

  combined_dat <- left_join(true_viral_loads_melted,obs_viral_loads_melted) %>% as_tibble

  return(combined_dat)
}

#'@export
simulate_infection_times <- function(n, prob_infection, overall_prob=NULL){
  if(is.null(overall_prob)){
    overall_prob <- sum(prob_infection)
  }
  scaled_prob<- prob_infection/sum(prob_infection)
  are_infected <- numeric(n)
  infection_times <- numeric(n)
  for(i in 1:n){
    infection <- rbinom(1,1, overall_prob)
    are_infected[i] <- infection
    if(infection == 1){
      t_inf <- sample(1:length(prob_infection), 1, prob=scaled_prob)
      infection_times[i] <- t_inf
    } else {
      infection_times[i] <- -1
    }
  }
  return(infection_times)
}

#' Simulate full line list data
#'
#' Takes a vector of infection incidence (absolute numbers) and returns a tibble with line list data for all
#' individuals in the population. Infected individuals are assigned symptom onset, an incubation periods (if applicable)
#' and confirmation delays from log normal and gamma distributions respectively.
#'
#' @export
simulate_observations_wrapper <- function(
  incidence, times, symp_frac=0.35,
  population_n=length(incidence),
  incu_period_par1=1.621,incu_period_par2=0.418,
  conf_delay_par1=5,conf_delay_par2=2){
  not_infected <- population_n-sum(incidence)
  inc_dat <- tibble(infection_time=c(NA,times),inc=c(not_infected,incidence))
  inc_dat <- inc_dat %>% uncount(inc)

  inc_dat <- inc_dat %>%
    mutate(i=1:n()) %>%
    mutate(is_infected = ifelse(is.na(infection_time), 0, 1)) %>%
    ## Is this individual going to be symptomatic?
    mutate(is_symp=ifelse(is_infected, rbinom(n(), 1, symp_frac), 0)) %>%
    ## Symptom onset time
    mutate(incu_period=ifelse(is_infected & is_symp, rlnorm(n(), incu_period_par1, incu_period_par2), NA),
           onset_time=infection_time+round(incu_period)) %>%
    ## Confirmation time
    mutate(confirmation_delay=extraDistr::rdgamma(n(),conf_delay_par1,conf_delay_par2))
  inc_dat
}
