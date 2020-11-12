#' @export
likelihood_cpp_wrapper <- function(obs_dat, ages, times,
                                          pars, prob_infection, pos_only=FALSE){
  liks_tj <- 0
  if(pos_only) {
    use_func <- likelihood_pos_only_cpp
  } else {
    use_func <- likelihood_cpp
  }
  for(i in seq_along(times)){
    ages1 <- ages[(times[i] - ages) > 0]
    obs1 <- obs_dat[[i]]
    liks_tj <- liks_tj + sum(use_func(obs1, times[i], ages1, pars, prob_infection))
  }
  liks_tj
}


#' Function to give probability of observing x given age a and the viral kinetics curve
#' @export
p_a <- function(x,a,pars,viral_loads) {
  viral_load_sd <- pars["obs_sd"]
  LOD <- pars["intercept"]
  renormalize <- extraDistr::pgumbel(LOD, viral_loads[a],viral_load_sd,lower.tail=TRUE)
  probs <- extraDistr::dgumbel(x, mu=viral_loads[a], sigma=viral_load_sd, log=FALSE)/renormalize
  probs
}

#' Probability of having a detectable Ct for a given time since infection
#' @export
prop_detectable <- function(a, pars,viral_loads){
  viral_load_sd <- pars["obs_sd"]
  LOD <- pars["intercept"]
  additional_prob <- 1

  ## How many days spent in detectability loss phase?
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  days_potential_loss <- a - t_switch
  if(days_potential_loss >= 0){
    additional_prob <- (1-pars["prob_detect"]*pars["t_unit"])^days_potential_loss
  }

  main_probs <- extraDistr::pgumbel(LOD,mu=viral_loads[a],sigma=viral_load_sd, lower.tail=TRUE, log.p=FALSE)
  main_probs * additional_prob
}

#'
#'@export
likelihood_R <- function(obs_dat, ages, pars, prob_infection){
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
    undetectable_prob <- 1 - sum(sapply(ages1, function(a) prob_detectable_dat[a]*prob_infection[obs_time-a]))

    detectable_prob <- function(X_i){
      sum(vapply(ages1,
                 FUN=function(a) p_a(X_i, a, pars, viral_loads)*prob_infection[obs_time-a]*prob_detectable_dat[a],
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
                 FUN=function(a) p_a(X_i, a, pars, viral_loads)*prob_infection[obs_time-a]*prob_detectable_dat[a],
                 FUN.VALUE=numeric(1)))
    }
    lik_tj <- lik_tj + log(sapply(obs1,function(x_i) detectable_prob(x_i)/prob_detectable_and_infected))
  }
  lik_tj
}


