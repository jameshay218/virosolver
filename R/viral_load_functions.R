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
  comb_dat <- NULL
  for(obs_time in obs_times){
    ages1 <- ages[(obs_time - ages) > 0]
    densities <- pred_dist_cpp(test_cts, ages1, obs_time, pars, prob_infection)
    comb_dat[[obs_time]] <- tibble(ct=test_cts,density=densities, t=obs_time)
  }
  comb_dat <- do.call("bind_rows",comb_dat)
  comb_dat

}
