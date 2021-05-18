#' @export
viral_load_func <- function(pars, obs_t, convert_vl=FALSE, infection_time=0){
  tshift <- pars["tshift"] + infection_time ## Days post infection until growth
  desired_mode <- pars["desired_mode"] ## Days post growth to peak
  t_switch <- pars["t_switch"] ## Days post peak to switch point

  height <- pars["viral_peak"] ## Peak Ct value
  level_switch <- pars["level_switch"] ## Ct value at the switch

  true_0 <- pars["true_0"] ## y-axis shift
  ct_intercept <- pars["intercept"] ## Ct intercept
  wane_rate <- (height-level_switch)/t_switch
  wane_rate2 <- (level_switch-ct_intercept)/pars["wane_rate2"]
  growth_rate <- (height-true_0)/desired_mode


  t_period1 <- obs_t <= tshift
  t_period2 <- obs_t > tshift & obs_t <= (desired_mode + tshift)
  t_period3 <- obs_t > (desired_mode + tshift) & obs_t <= (desired_mode + tshift + t_switch)
  t_period4 <- obs_t > (desired_mode + tshift + t_switch)
  y <- numeric(length(obs_t))
  y[t_period1] <- true_0
  y[t_period2] <- growth_rate * (obs_t[t_period2] - tshift) + true_0
  y[t_period3] <- height - wane_rate * (obs_t[t_period3] - (desired_mode+tshift))
  y[t_period4] <- level_switch - wane_rate2*(obs_t[t_period4] - (desired_mode + tshift + t_switch))
  ct <- y
  if(convert_vl){
    vl <- ((ct_intercept-ct)/log2(10)) + pars["LOD"]
    y <- vl
  }
  y
}

#' @export
pred_dist_wrapper <- function(test_cts, obs_times, ages, pars, prob_infection){
  ## Because we have a different standard deviation for different
  ## Time at which standard deviation is reduced
  t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
  sd_mod <- rep(pars["sd_mod"], length(ages)) ## Up until t_switch, full standard deviation

  ## Prior to t_switch, 1
  unmod_vec <- 1:min(t_switch,length(ages))
  sd_mod[unmod_vec] <- 1

  ## For the next sd_mod_wane days, variance about modal Ct trajectory decrease linearly
  decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
  ## The rest are at sd_mod

  comb_dat <- NULL
  for(obs_time in obs_times){
    ages1 <- ages[(obs_time - ages) > 0]
    densities <- pred_dist_cpp(test_cts, ages1, obs_time, pars, prob_infection,sd_mod)
    comb_dat[[obs_time]] <- tibble(ct=test_cts,density=densities, t=obs_time)
  }
  comb_dat <- do.call("bind_rows",comb_dat)
  comb_dat

}


#' @export
simulate_viral_loads_example <- function(ages, kinetics_pars,N=100){
  t_switch <-  kinetics_pars["t_switch"] + kinetics_pars["desired_mode"] + kinetics_pars["tshift"]
  sd_mod <- rep(kinetics_pars["sd_mod"], max(ages))
  unmod_vec <- 1:min(t_switch,max(ages))
  sd_mod[unmod_vec] <- 1
  decrease_vec <- (t_switch+1):(t_switch+kinetics_pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-kinetics_pars["sd_mod"])/kinetics_pars["sd_mod_wane"])*seq_len(kinetics_pars["sd_mod_wane"])
  
  sim_dat <- matrix(ncol=N,nrow=length(ages))
  for(age in ages){
    ## For each value we're going to simulate, pre-determine if it will still be detectable or not
    detectable_statuses <- rnbinom(N, 1, prob=kinetics_pars["prob_detect"]) + 
      kinetics_pars["tshift"] + kinetics_pars["desired_mode"] + kinetics_pars["t_switch"]
    cts <- rep(virosolver::viral_load_func(kinetics_pars, age, FALSE, 0),N)
    cts[detectable_statuses <= age] <- 1000
    sd_used <- kinetics_pars["obs_sd"]*sd_mod[age]
    ct_obs_sim <- extraDistr::rgumbel(N, cts, sd_used)
    ct_obs_sim <- pmin(ct_obs_sim, kinetics_pars["intercept"])
    sim_dat[age,] <- ct_obs_sim
  }
  sim_dat <- data.frame(sim_dat)
  colnames(sim_dat) <- 1:N
  sim_dat$time_since_infection <- ages
  sim_dat <- sim_dat %>% pivot_longer(-time_since_infection) %>% rename(ct=value,age=time_since_infection,i=name)
  return(sim_dat)
}
