#' @export
viral_load_func <- function(pars, obs_t, convert_ct=TRUE, infection_time=0){
  tshift <- pars["tshift"] + infection_time ## Days post infection until growth
  desired_mode <- pars["desired_mode"] #+ infection_time## Days post growth to peak
  t_switch <- pars["t_switch"] ## Days post peak to switch point

  height <- pars["viral_peak"] ## Peak viral load
  obs_sd <- pars["obs_sd"] ## Variance about mean
  level_switch <- pars["level_switch"]

  true_0 <- pars["true_0"] ## y-axis shift
  yintercept <- pars["intercept"] ## Ct intercept

  wane_rate <- (height-level_switch)/t_switch
  wane_rate2 <- (level_switch - pars["LOD"])/pars["wane_rate2"]
  print(pars["wane_rate2"])
  growth_rate <- (height-true_0)/desired_mode

  y <- numeric(length(obs_t))
  y[obs_t <= tshift] <- true_0
  y[obs_t > tshift & obs_t <= (desired_mode + tshift)] <- growth_rate * (obs_t[obs_t > tshift & obs_t <= (desired_mode + tshift)] - tshift) + true_0
  y[obs_t > (desired_mode + tshift) & obs_t <= (desired_mode + tshift + t_switch)] <- height - wane_rate * (obs_t[obs_t > (desired_mode + tshift) & obs_t <= (desired_mode + tshift + t_switch)] - (desired_mode+tshift))
  y[obs_t > (desired_mode + tshift + t_switch)] <- level_switch - wane_rate2*(obs_t[obs_t > (desired_mode + tshift + t_switch)] - (desired_mode + tshift + t_switch))
  ct <- y
  if(convert_ct){
    ct <- yintercept - log2(10) * (ct-pars["LOD"])
  }
  ct
}
