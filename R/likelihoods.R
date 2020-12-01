#' @export
likelihood_cpp_wrapper <- function(obs_dat, ages, times,
                                   pars, prob_infection,
                                   pos_only=FALSE,
                                   undetectable_counts=NULL){
  ## If flag to only use detectable Ct values or not, use different likelihood function
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
  ## Max sure we don't go past maximum of ages vector
  unmod_vec <- 1:min(t_switch,max(ages)+1)
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
  ## The rest are at sd_mod

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

    liks_tj <- liks_tj +
      sum(use_func(obs1, times[i], ages1, pars, prob_infection,sd_mod)) + ## Detectable Cts
      undetectable_lik
  }
  liks_tj
}


#' Function to give probability of observing x given age a and the viral kinetics curve
#' @export
p_a <- function(x,a,pars,viral_loads,sd_mod) {
  viral_load_sd <- pars["obs_sd"]*sd_mod[a]
  LOD <- pars["intercept"]
  renormalize <- extraDistr::pgumbel(LOD, viral_loads[a],viral_load_sd,lower.tail=TRUE) -
    extraDistr::pgumbel(0, viral_loads[a],viral_load_sd,lower.tail=TRUE)
  probs <- extraDistr::dgumbel(x, mu=viral_loads[a], sigma=viral_load_sd, log=FALSE)/renormalize
  probs
}

#' Probability of having a detectable Ct for a given time since infection
#' @export
prop_detectable_single <- function(a, pars,viral_loads, sd_mod){
  viral_load_sd <- pars["obs_sd"]*sd_mod[a]
  LOD <- pars["intercept"]
  additional_prob <- 1
  ## How many days spent in detectability loss phase?
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]

  days_potential_loss <- a - t_switch
  if(days_potential_loss >= 0){
    additional_prob <- (1-pars["prob_detect"]*pars["t_unit"])^days_potential_loss
  }
  main_probs <- extraDistr::pgumbel(LOD,mu=viral_loads[a],sigma=viral_load_sd, lower.tail=TRUE, log.p=FALSE)
  #main_probs <- pnorm(LOD,mean=viral_loads[a],sd=viral_load_sd, lower.tail=TRUE, log.p=FALSE)
  main_probs * additional_prob
}

#' @export
prop_detectable <- function(a, pars,viral_loads){
  ## Because we have a different standard deviation for different
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], max(ages))

  ## Prior to t_switch, full standard deviation
  ## Max sure we don't go past maximum of ages vector
  unmod_vec <- 1:min(t_switch,max(ages)+1)
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
  ## The rest are at sd_mod

  vapply(a, function(x){
    prop_detectable_single(x, pars, viral_loads, sd_mod)
  },FUN.VALUE=numeric(1))
}

#'
#'@export
likelihood_R <- function(obs_dat, ages, pars, prob_infection){
  LOD <- pars["intercept"]

  ## Because we have a different standard deviation for different
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], length(ages)) ## Up until t_switch, full standard deviation

  ## Prior to t_switch, 1
  unmod_vec <- 1:min(t_switch,length(ages))
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
  ## The rest are at sd_mod

  viral_loads <- viral_load_func(pars, ages)
  prob_detectable_dat <- sapply(ages, function(a) prop_detectable_cpp(a, viral_loads[a],
                                                                      pars["obs_sd"]*sd_mod[a], pars["intercept"],
                                                                      t_switch, pars["prob_detect"]))

  times <- unique(obs_dat$t)
  lik_tj <- 0
  for(obs_time in times){
    ages1 <- ages[(obs_time - ages) > 0]
    obs1 <- obs_dat %>% filter(t == obs_time) %>% pull(ct)
    undetectable_prob <- 1 - sum(sapply(ages1, function(a) prob_detectable_dat[a]*prob_infection[obs_time-a]))

    detectable_prob <- function(X_i){
      sum(vapply(ages1,
                 FUN=function(a) p_a(X_i, a, pars, viral_loads,sd_mod)*prob_infection[obs_time-a]*prob_detectable_dat[a],
                 FUN.VALUE=numeric(1)))
    }
    lik_tj <- lik_tj + log(sapply(obs1,function(x_i) (x_i >= LOD)*undetectable_prob + (x_i < LOD)*detectable_prob(x_i)))
  }
  lik_tj
}


#'@export
likelihood_pos_R <- function(obs_dat, ages, pars, prob_infection){
  LOD <- pars["intercept"]
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]

  viral_loads <- viral_load_func(pars, ages)
  prob_detectable_dat <- sapply(ages, function(a) prop_detectable_cpp(a, viral_loads[a],
                                                                      pars["obs_sd"], pars["intercept"],
                                                                      t_switch, pars["prob_detect"]))
  times <- unique(obs_dat$t)
  lik_tj <- 0
  for(obs_time in times){
    ages1 <- ages[(obs_time - ages) > 0]
    obs1 <- obs_dat %>% filter(t == obs_time) %>% pull(ct)
    prob_detectable_and_infected <- sum(sapply(ages1, function(a) prob_detectable_dat[a]*prob_infection[obs_time-a]))

    detectable_prob <- function(X_i){
      sum(vapply(ages1,
                 FUN=function(a) p_a(X_i, a, pars, viral_loads,sd_mod)*prob_infection[obs_time-a]*prob_detectable_dat[a],
                 FUN.VALUE=numeric(1)))
    }
    lik_tj <- lik_tj + log(sapply(obs1,function(x_i) detectable_prob(x_i)/prob_detectable_and_infected))
  }
  lik_tj
}


