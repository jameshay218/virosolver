#' Calculate GR CT relationship for one R0
#' @export
#' @examples
#' times <- 0:150 # Duration of epidemic in days
#' lastday <- 35 # For each individual, how long do we assume their measurements can contribute towards the model? 
#' cts <- seq(0,40,by=1)
#' model_pars <- c("R0"=R0,"infectious"=5,"incubation"=2,"I0"=0.0001,"t0"=0)
#' viral_load_pars <- c(beta = 0.1, tshift = 0, desired_mode = 5, viral_peak = 19.7, 
#'                     obs_sd = 5, sd_mod = 0.789, sd_mod_wane = 14, true_0 = 40, 
#'                     intercept = 40, LOD = 3, incu = 5, t_switch = 13.3, level_switch = 38, 
#'                     wane_rate2 = 1000, prob_detect = 0.103, t_unit = 1)
#' calculate_gr_ct_relationships(model_pars, viral_load_pars, times, lastday, cts)
calculate_gr_ct_relationships <- function(model_pars, viral_load_pars, times, lastday=35,cts=seq(0,40,by=1),population_n=1,convert_daily=TRUE,tstep=1){
  skew <- function(values,weights) {
    weights_std <- weights/sum(weights, na.rm=TRUE)
    xbar <- sum(values*weights_std, na.rm=TRUE)
    xi_xbar <- values - xbar
    return((sum(weights_std*xi_xbar^3))/((sum(weights_std*xi_xbar^2))^(3/2)))
  }
  epidemic_process <- virosolver::simulate_seir_wrapper(population_n, times, model_pars,convert_daily = convert_daily)
  if(convert_daily) times <- unique(floor(times))
  #times <- unique(floor(times))
  incidence <- epidemic_process$per_cap_incidence
  gr <- epidemic_process$growth_rates
  #### Calculating viral load and age distributions: ####
  age_res_detectable <- matrix(0,nrow=length(times),ncol=lastday)
  age_res <- matrix(0,nrow=length(times),ncol=lastday)
  
  # Detectable probabilities over [lastday] days prior to test
  viral_loads <- viral_load_func(viral_load_pars, 1:lastday)
  detectable_props <- prop_detectable(1:lastday, viral_load_pars, viral_loads)
  ages <- 1:lastday

  # Create matrices of relative frequencies of observing age of infection and viral load by day of testing:
  dat <- pred_dist_wrapper(cts, seq_along(times),1:lastday,viral_load_pars, incidence)
  
  for (i in 2:length(times)) {
    past_inc <- incidence[(i-1):(max(i-lastday,1))]
    days <- 1:length(past_inc)
    age_res_detectable[i,days] <- past_inc*detectable_props[days]
    age_res[i,days] <- past_inc
  }

  # Summaries of ages from detectable individuals
  age_res_std <- age_res_detectable/apply(age_res_detectable, 1, sum, na.rm=TRUE)
  age_mean <- apply(age_res_std, 1, function(res) sum(res*(1:lastday), na.rm=TRUE))
  age_skew <- apply(age_res_std, 1, function(row) skew(values=(1:lastday), weights=row))
  age_res_std_csum <- t(apply(age_res_std, 1, FUN=function(res) cumsum(res)))
  age_median <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.5]))
  age_lower <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.25]))
  age_upper <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.75]))
  
  age_data <- tibble(t=times, age_mean=age_mean, age_median=age_median, 
                     age_skew=age_skew,age_lower25=age_lower,age_upper75=age_upper)
  
  # Summaries of ages from all individuals
  age_res_std <- age_res/apply(age_res, 1, sum, na.rm=TRUE)
  age_mean <- apply(age_res_std, 1, function(res) sum(res*(1:lastday), na.rm=TRUE))
  age_skew <- apply(age_res_std, 1, function(row) skew(values=(1:lastday), weights=row))
  age_res_std_csum <- t(apply(age_res_std, 1, FUN=function(res) cumsum(res)))
  age_median <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.5]))
  age_lower <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.25]))
  age_upper <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.75]))
  
  age_data_all <- tibble(t=times, age_mean=age_mean, age_median=age_median, 
                     age_skew=age_skew,age_lower25=age_lower,age_upper75=age_upper)
  
  p_ages <- ggplot(age_data) +
    geom_line(aes(x=t,y=age_mean,col="Mean")) +
    geom_line(aes(x=t,y=age_median,col="Median")) +
    ylab("Time since infection") +
    xlab("Time since start of epidemic") +
    scale_color_manual(name="Metric",values=c("Mean"="purple","Median"="orange")) +
    theme_bw()
  
  p_ages_all <- ggplot(age_data_all) +
    geom_line(aes(x=t,y=age_mean,col="Mean")) +
    geom_line(aes(x=t,y=age_median,col="Median")) +
    ylab("Time since infection") +
    xlab("Time since start of epidemic") +
    scale_color_manual(name="Metric",values=c("Mean"="purple","Median"="orange")) +
    theme_bw()
  
  
  p_ages_skew <- ggplot(age_data) +
    geom_line(aes(x=t,y=age_skew,col="Skew")) +
    ylab("Time since infection (skew)") +
    xlab("Time since start of epidemic") +
    scale_color_manual(name="Metric",values=c("Mean"="purple","Median"="orange","Skew"="darkgreen")) +
    theme_bw()
  
  ## Get mean of detectable Ct values
  ct_dat_mean <- dat %>% filter(ct < max(cts)) %>% group_by(t) %>% 
    mutate(density_scaled=density/sum(density)) %>% 
    summarize(mean_ct=sum(ct*density_scaled)) 
  
  ## Get median of detectable Ct values
  ct_dat_median <- dat %>% filter(ct < max(cts)) %>% group_by(t) %>% 
    mutate(density_scaled=density/sum(density)) %>% 
    mutate(cumu_density=cumsum(density_scaled)) %>%
    filter(cumu_density >= 0.5) %>%
    filter(ct == min(ct)) %>% 
    select(ct, t) %>%
    rename(median_ct = ct)
  
  ct_dat_lower25 <- dat %>% filter(ct < max(cts)) %>% group_by(t) %>% 
    mutate(density_scaled=density/sum(density)) %>% 
    mutate(cumu_density=cumsum(density_scaled)) %>%
    filter(cumu_density >= 0.25) %>%
    filter(ct == min(ct)) %>% 
    select(ct, t) %>%
    rename(median_ct = ct)
  
  ct_dat_upper75 <- dat %>% filter(ct < max(cts)) %>% group_by(t) %>% 
    mutate(density_scaled=density/sum(density)) %>% 
    mutate(cumu_density=cumsum(density_scaled)) %>%
    filter(cumu_density >= 0.75) %>%
    filter(ct == min(ct)) %>% 
    select(ct, t) %>%
    rename(median_ct = ct)

  ## Get skew of detectable Ct values
  ct_skew <- dat %>% filter(ct < max(cts)) %>%
    group_by(t) %>% 
    mutate(density_scaled=density/sum(density)) %>% 
    mutate(mean_ct=sum(ct*density_scaled),top_part=(ct-mean_ct)^3,bot_part=(ct-mean_ct)^2) %>% 
    group_by(t) %>% 
    summarize(top=sum(top_part*density_scaled),bot=(sum(bot_part*density_scaled)^1.5),skew=top/bot) %>%
    select(t, skew) %>%
    rename(skew_ct = skew)
  
  ct_dat <- left_join(ct_dat_mean, ct_dat_median, by="t") %>% 
    left_join(ct_skew, by="t") %>%
    left_join(ct_dat_lower25, by="t") %>% 
    left_join(ct_dat_upper75, by="t")
  
  ## Plot mean of detectable Cts
  p_ct <- ggplot(ct_dat) + 
    geom_line(aes(x=t,y=mean_ct,col="Mean")) +
    geom_line(aes(x=t,y=median_ct, col="Median")) +
    scale_color_manual(name="Metric",values=c("Mean"="red","Median"="blue")) +
    theme_bw() +
    ylab("Average Ct value") +
    xlab("Time since start of epidemic")
  
  all_dat <- left_join(gr, age_data, by="t") %>% left_join(ct_dat, by="t") %>% 
    left_join(epidemic_process$seir_outputs %>% rename(t=step), by="t") %>%
    left_join(tibble(t=unique(floor(times)),inc=epidemic_process$per_cap_inc), by="t")
  
  all_dat_final <- all_dat %>% dplyr::filter(t>min(times),t<max(times)) %>% 
    select(ver, t, GR, age_mean,age_median,age_skew,mean_ct,median_ct,skew_ct,Rt,inc)
  
  all_dat <- all_dat %>% dplyr::filter(ver == "daily",t>min(times),t<max(times)) %>% 
    select(t, GR, age_mean,age_median,age_skew,mean_ct,median_ct,skew_ct,Rt,inc)
  
  theme_use <- theme_bw()
  scale_color_use <- scale_color_viridis_c()
  max_GR <- max(abs(all_dat$GR),na.rm=TRUE)
  x_axis_use <- scale_x_continuous(limits=c(-max_GR,max_GR))
  
  p1 <- ggplot(all_dat) +
    geom_path(aes(x=GR,y=age_mean,col=t),linewidth=0.75) +
    geom_vline(xintercept=0,linetype="dashed") +
    theme_use+
    scale_color_use +
    xlab("Growth rate") +
    ylab("Mean time-since-infection")+
    x_axis_use
  p2 <-  ggplot(all_dat) +
    geom_path(aes(x=GR,y=age_skew,col=t),linewidth=0.75) +
    geom_vline(xintercept=0,linetype="dashed") +
    xlab("Growth rate") +
    ylab("Skew of time-since-infection")+
    theme_use +
    scale_color_use+
    x_axis_use
  p3 <-  ggplot(all_dat) +
    geom_path(aes(x=GR,y=mean_ct,col=t),linewidth=0.75) +
    geom_vline(xintercept=0,linetype="dashed") +
    xlab("Growth rate") +
    ylab("Mean Ct value")+
    scale_y_continuous(trans="reverse")+
    theme_use+
    scale_color_use+
    x_axis_use
  
  p4 <-  ggplot(all_dat) +
    geom_path(aes(x=GR,y=skew_ct,col=t),linewidth=0.75) +
    geom_vline(xintercept=0,linetype="dashed") +
    xlab("Growth rate") +
    ylab("Skew of Ct values")+
    theme_use+
    scale_color_use+
    x_axis_use
  
  
  
  p_main <- p1 + p2 + p3 + p4 + plot_layout(guides="collect")
  
  
  return(list(p_ct=p_ct, 
              p_ages_detectable=p_ages,
              p_ages_all=p_ages_all,
              growth_rates=gr,
             tsi_detectable=age_data, 
              tsi_all=age_data_all,
              ct_values=dat,
             p_main=p_main,
             all_dat=all_dat_final))
}

#' Calculate GR CT relationship for multiple R0s
#' @export
calculate_gr_ct_relationships_R0s <- function(model_pars, viral_load_pars, times, lastday=35, cts=seq(0,40,by=0.1), R0s){
  ## Solve model for each R0
  for(R0 in R0s){
    model_pars["R0"] <- R0
    res <- virosolver::calculate_gr_ct_relationships(model_pars,viral_load_pars,times,lastday,cts,convert_daily=TRUE)
    res$all_dat$R0 <- R0
    all_dat <- bind_rows(all_dat, res$all_dat)
  }
  
  all_dat <- all_dat %>% filter(ver == "daily")
  
  ## Age plots
  p_rt_age_skew <- ggplot(all_dat %>% 
                            filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=Rt,y=age_skew,col=R0,group=R0)) +
    geom_vline(xintercept=1,linetype="dashed") +
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Skew of TSI") +
    xlab("Rt") +
    theme_use
  p_rt_age_mean <- ggplot(all_dat %>% drop_na() %>% 
                            filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=Rt,y=age_mean,col=R0,group=R0)) +
    geom_vline(xintercept=1,linetype="dashed") +
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Mean of TSI") +
    xlab("Rt") +
    theme_use
  
  p_gr_age_skew <- ggplot(all_dat %>% 
                            filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=GR,y=age_skew,col=R0,group=R0)) +
    geom_vline(xintercept=0,linetype="dashed") +
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Skew of TSI") +
    xlab("Growth rate of infections") +
    theme_use
  
  p_gr_age_mean <- ggplot(all_dat %>% drop_na() %>% 
                            filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=GR,y=age_mean,col=R0,group=R0)) +
    geom_vline(xintercept=0,linetype="dashed") +
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Mean of TSI") +
    xlab("Growth rate of infections") +
    theme_use
  
  p_age <- p_rt_age_mean + p_rt_age_skew + p_gr_age_mean + p_gr_age_skew + plot_layout(guides="collect")
  
  ## Ct plots
  p_rt_mean <- ggplot(all_dat %>% drop_na() %>% 
                        filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=Rt,y=mean_ct,col=R0,group=R0)) +
    geom_vline(xintercept=1,linetype="dashed") +
    scale_color_gradient(low="blue",high="red")+
    ylab("Mean of Ct values") +
    xlab("Rt") +
    theme_use
  
  p_rt_skew <- ggplot(all_dat %>% filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=Rt,y=skew_ct,col=R0,group=R0)) +
    geom_vline(xintercept=1,linetype="dashed") +
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Skew of Ct values") +
    xlab("Rt") +
    theme_use
  
  p_gr_mean <- ggplot(all_dat %>% drop_na()%>% 
                        filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=GR,y=mean_ct,col=R0,group=R0)) +
    geom_vline(xintercept=0,linetype="dashed") +
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Mean of Ct values") +
    xlab("Growth rate of infections") +
    theme_use
  
  p_gr_skew <- ggplot(all_dat %>% 
                        filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=GR,y=skew_ct,col=R0,group=R0)) +
    geom_vline(xintercept=0,linetype="dashed") +
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    scale_x_continuous(trans="reverse") +
    ylab("Skew of Ct values") +
    xlab("Growth rate of infections") +
    theme_use
  
  p_cts <- p_rt_mean + p_rt_skew + p_gr_mean + p_gr_skew + plot_layout(guides="collect")
  
  
  p_both_ctA <- ggplot(all_dat %>% 
                         filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=mean_ct,y=skew_ct,col=R0,group=R0))+
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Skew of Ct values") +
    xlab("Mean of Ct values") +
    theme_use 
  
  p_both_age <- ggplot(all_dat %>% 
                         filter(inc > model_pars["I0"])) + 
    geom_path(aes(x=age_mean,y=age_skew,col=R0,group=R0))+
    scale_color_gradient(low="blue",high="red",limits=range(R0s))+
    ylab("Skew of TSI") +
    xlab("Mean of TSI") +
    theme_use 
  
  p_tsi_dist <- ggplot(all_dat %>% filter(t < lastday) %>% group_by(R0) %>% 
                         mutate(inc=inc/sum(inc))) + 
    geom_line(aes(x=lastday-t,y=inc,group=R0,col=R0)) +
    geom_vline(data=all_dat %>% filter(t == lastday),aes(xintercept=age_median)) +
    scale_color_gradient(low="blue",high="red",limits=range(R0s)) + 
    xlab("Time since infection") +
    ylab("Density") +
    ggtitle(paste0("Time-since-infection distribution at day ", lastday))+
    facet_wrap(~R0)
  
  return(list(p_age, p_cts, p_both_ctA, p_both_age, p_tsi_dist, all_dat))
}

#' Compare TSI vs. GR
#' @export
compare_tsi_dist_at_gr <- function(dat, grs,tolerance=0.01,lastday=35){
  
  get_tsi_dist_at_gr <- function(dat, gr=0){
    dat %>% 
      filter(inc > sir_pars["I0"]) %>% 
      ungroup() %>% 
      mutate(diff_gr=abs(GR - gr)) %>% 
      group_by(R0) %>% 
      mutate(use=diff_gr==min(diff_gr)) %>%
      filter(use == TRUE) %>%
      select(t, R0, GR, diff_gr) %>%
      ungroup()
  }
  use_tsi <- NULL
  for(gr in grs){
    tmp <- get_tsi_dist_at_gr(dat, gr) %>% 
      rename(use_t = t, closest_GR=GR) %>% 
      filter(diff_gr < tolerance) %>%
      mutate(GR_use=gr)
    use_tsi <- bind_rows(use_tsi, left_join(dat %>% ungroup(),tmp,by="R0") %>% 
                           filter(!is.na(use_t), t <= use_t, t > use_t - lastday) %>%
                           group_by(R0) %>% 
                           mutate(t = t-max(t)) %>%
                           group_by(R0) %>% 
                           mutate(inc = inc/sum(inc)))
  }
  
  
  p1 <- ggplot(use_tsi) + 
    geom_line(aes(x=-t, y=inc,col=R0,group=R0),linewidth=0.75) +
    scale_color_viridis_c() +
    xlab("Time since infection") +
    ylab("Density") +
    facet_wrap(~GR_use,scales="free_y")
  p1
}
