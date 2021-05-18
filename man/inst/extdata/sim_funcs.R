## Simulate deterministic SEIR model
## INPUTS: Takes a vector of SEIR model parameters, times to solve over and population size.
##      4. switch_model: if TRUE, uses the SEIR model with 2 switch points in transmission intensity
## OUTPUTS: 
##      1. Plot of all SEIR compartments over time
##      2. Plot of incidence and prevalence over time
##      3. Absolute incidence per time point
##      4. Per capita incidence per time point
##      5. Per capita prevalence (E+I) per time point
##      6. Overall probability of infection
simulate_seir_process <- function(pars, times, N=1){
  ## Pull parameters for SEIR model
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
  ## Set up initial conditions.
  ## Note if population_n=1, then solves everything per capita
  # init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0)
  init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0,0)
  
  ## Solve the SEIR model using the rlsoda package
  sol <- rlsoda::rlsoda(init, times, C_SEIR_model_rlsoda, parms=seir_pars, dllname="virosolver",
                        deSolve_compatible = TRUE,return_time=TRUE,return_initial=TRUE,atol=1e-10,rtol=1e-10)
    

  ## Convert to data frame and column names
  sol <- as.data.frame(sol)
  colnames(sol) <- c("time","S","E","I","R","cumu_exposed","cumu_infected")
  ## Get Rt
  sol$Rt <- (sol$S) * pars["R0"]
  
  ## Shift for start
  sol$time <- sol$time + floor(pars["t0"])
  
  ## Dummy rows from pre-seeding
  if(pars["t0"] > 0){
    dummy_row <- data.frame("time"=0:(floor(unname(pars["t0"]))-1),"S"=N,"E"=0,"I"=0,"R"=0,"cumu_exposed"=0,"cumu_infected"=0,"Rt"=unname(pars["R0"])*N)
    sol <- bind_rows(dummy_row, sol)
  }
  sol <- sol[sol$time %in% times,]
  
  ## Pull raw incidence in absolute numbers
  inc <- c(0,diff(sol[,"cumu_exposed"]))
  inc <- pmax(0, inc)
  
  ## Get per capita incidence and prevalence
  per_cap_inc <- inc/N
  per_cap_prev <- (sol$E + sol$I)/N
  
  ## Melt solution, get per capita and plot
  sol <- reshape2::melt(sol, id.vars="time")
  sol$value <- sol$value/N
  
  ## Plot all compartments
  p <- ggplot(sol) +
    geom_line(aes(x=time,y=value,col=variable)) +
    ylab("Per capita") +
    xlab("Date") +
    theme_bw()
  
  ## Plot incidence and prevalence
  p_inc <- ggplot(data.frame(x=times,y=per_cap_inc,y1=per_cap_prev)) +
    geom_line(aes(x=x,y=y),col="red") +
    geom_line(aes(x=x,y=y1),col="blue") +
    ylab("Per capita incidence (red) and prevalence (blue)") +
    xlab("Date") +
    theme_bw()
  
  return(list(plot=p,
              incidence_plot=p_inc,
              seir_outputs=sol,
              raw_incidence=inc,
              per_cap_incidence=per_cap_inc,
              per_cap_prevalence=per_cap_prev,
              overall_prob_infection=sum(per_cap_inc)))
}


## Subset line list data by testing strategy. Options:
#' 1. Sample a random fraction of the population if the only argument is frac_report
#' 2. Sample some random fraction of the population at a subset of time points, specified by timevarying_prob
#' 3. Observe symptomatic individuals with some fixed probability, frac_report if symptomatic is TRUE
#' 4. Observe symptomatic individuals with some time-varying probability, timevarying_prob, if symptomatic is TRUE
#' INPUTS: 
#'      1. individuals: the full line list from the simulation, returned by virosolver::simulate_observations_wrapper
#'      2. solve_times: vector of times at which individuals can be reported
#'      3. frac_report: the overall fraction/probability of individuals who are reported
#'      4. timevarying_prob: a tibble with variables t and prob. This gives the probability of being reported on day t
#'      5. symptomatic: if TRUE, then individuals are reported after developing symptoms. If FALSE, then we take a random cross-section 
#' OUTPUTS: 
#'      1. A tibble with line list data for individuals who were observed
#'      2. A plot of incidence for both observed individuals and the entire simulated population
#'      3. Plot growth rate of cases/infections in the entire population and observed population
simulate_reporting <- function(individuals,
                               solve_times, 
                               frac_report=1,
                               timevarying_prob=NULL, 
                               symptomatic=FALSE){
  ## The basic form is that a constant fraction of individuals are reported each day
  n_indivs <- nrow(individuals)
  sampled_individuals <- individuals
  
  ## If a flat reporting rate
  if(is.null(timevarying_prob)){
    ## Base case, just sampled individuals at random
    if(!symptomatic){
      ## Sample a fraction of the overall line list and assign each individual a sample time (at random)
      ## and a confirmation time (with confirmation delay)
      sampled_individuals <- sampled_individuals %>% 
        sample_frac(frac_report) %>%
        group_by(i) %>%
        mutate(sampled_time=sample(solve_times,n()),
               ## We sample but then have to wait for a result
               confirmed_time=sampled_time+confirmation_delay) %>%
        ungroup()
    } else {
      ## Symptomatic based surveillance. Observe individuals based on symptom onset date
      ## Subset by symptomatic individuals, observe some fraction of these at some delay post symptom onset
      sampled_individuals <- sampled_individuals %>% 
        filter(is_symp==1) %>% 
        sample_frac(frac_report) %>%
        mutate(sampled_time=onset_time+confirmation_delay,
               ## We sample individuals some number of days after symptom onset
               confirmed_time=sampled_time)
    }
    ## Reporting rate varies over time
  } else {
    ## Sample entire population
    if(!symptomatic){
      indivs <- unique(sampled_individuals$i)
      indivs_test <- indivs
      tmp_sampled_indivs <- NULL
      for(index in 1:nrow(timevarying_prob)) {
        sample_n <- round(timevarying_prob$prob[index]*length(indivs))
        sampled_time1 <- timevarying_prob$t[index]
        sampled_indivs <- sample(indivs_test, sample_n,replace=FALSE)
        tmp_sampled_indivs[[index]] <- sampled_individuals %>% 
          filter(i %in% sampled_indivs) %>% 
          mutate(sampled_time=sampled_time1,
                 confirmed_time = sampled_time + confirmation_delay)
        indivs_test <- setdiff(indivs_test, sampled_indivs)
      }
      sampled_individuals <- do.call("bind_rows", tmp_sampled_indivs)
      
      #browser()
      ## How many individuals are we going to observe by the end?
      #frac_report_overall <- 1-prod(1-timevarying_prob$prob)
      #frac_report_overall <- sum(timevarying_prob$prob)
      ## We will first get the fraction of individuals we will sample over the whole period
      ## Then, we will change the time-varying reporting probability to the relative number
      ## sampled on each day
      #scaled_timevarying_prob <- timevarying_prob$prob/frac_report_overall
      
      ## On each time point in timevarying_prob$t, choose some random fraction of the population
      ## to observe, weighted by scaled_timevarying_prob
      ## assign a sample time and confirmation time
      #sampled_individuals <- sampled_individuals %>%
      #  sample_frac(frac_report_overall) %>%
      # group_by(i) %>%
      #  mutate(sampled_time = sample(timevarying_prob$t, n(), replace=TRUE, prob=scaled_timevarying_prob),
      #        confirmed_time = sampled_time + confirmation_delay) %>%
      #  ungroup()
    } else {
      ## This is quite different - if you have symptom onset on day t, there is a probability that you will be observed
      ## Symptomatic based surveillance. Observe individuals based on symptom onset date
      sampled_individuals <- sampled_individuals %>% 
        filter(is_symp==1) %>% ## Subset for symptomatic
        left_join(timevarying_prob %>% rename(onset_time=t) %>% ## Join with time-varying reporting rate table
                    dplyr::select(-ver), by="onset_time") %>%
        group_by(i) 
      
      ## Quicker to vectorize
      sampled_individuals$is_reported <- rbinom(nrow(sampled_individuals), 1, sampled_individuals$prob) ## Assign individuals as detected or not
      
      sampled_individuals <- sampled_individuals %>%
        filter(is_reported == 1) %>% ## Only take detected individuals
        mutate(sampled_time=onset_time+confirmation_delay,
               ## We sample individuals some number of days after symptom onset
               confirmed_time=sampled_time)
    }
  }
  
  ## Plot incidence of infections,  onsets, confirmations and number sampled per day
  ## Get grouped (not line list) subset data
  grouped_dat <- sampled_individuals %>% 
    dplyr::select(i, infection_time, onset_time, sampled_time, confirmed_time) %>%
    pivot_longer(-i) %>% 
    drop_na() %>%
    group_by(name, value) %>% 
    tally() %>%
    rename(var=name,
           t=value,
           n=n) %>%
    complete(var, nesting(t),fill = list(n = 0)) %>%
    mutate(ver="Sampled individuals")
  
  ## Grouped entire dataset
  grouped_dat_all <- individuals %>% 
    dplyr::select(i, infection_time, onset_time) %>%
    pivot_longer(-i) %>% 
    drop_na() %>%
    group_by(name, value) %>% 
    tally() %>%
    rename(var=name,
           t=value,
           n=n) %>%
    complete(var, nesting(t),fill = list(n = 0)) %>%
    mutate(ver="All individuals")
  
  grouped_dat_combined <- bind_rows(grouped_dat, grouped_dat_all)
  
  ## Plot line list data
  p_all <- grouped_dat_combined %>% ggplot() +
    geom_line(aes(x=t,y=n,col=var)) +
    theme_bw() +
    ylab("Number of individuals") +
    xlab("Time") +
    theme(legend.position="bottom") +
    facet_wrap(~ver,ncol=1, scales="free_y")
  
  ## Get day-by-day growth rate
  grouped_dat_combined <- grouped_dat_combined %>% group_by(var, ver) %>% 
    mutate(gr_daily=log(n/lag(n,1))) %>%
    mutate(gr_daily = ifelse(is.finite(gr_daily), gr_daily, NA)) %>%
    ungroup()
  
  ## Get rolling average growth rate
  growth_rates_all <- expand_grid(grouped_dat_combined, window=seq(10,50,by=10)) %>% 
    arrange(var, ver, window, t) %>%
    group_by(var, ver, window) %>%
    mutate(gr_window=zoo::rollmean(gr_daily, window,align="right",fill=NA)) %>%
    mutate(window = as.factor(window))
  
  p_gr <- ggplot(growth_rates_all) +
    geom_line(aes(x=t,y=gr_window,col=ver)) +
    geom_hline(yintercept=0,linetype="dashed") +
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    xlab("Date") +    
    theme_bw() +
    theme(legend.position="bottom") +
    facet_grid(window~var)
  
  
  list(sampled_individuals=sampled_individuals,
       plot=p_all, 
       plot_gr=p_gr)
}

## Simulate observed Ct values for the line list dataset. NOTE this differs to virosolver::simulate_viral_loads,
## as this function only solves the viral kinetics model for the observation time
#' INPUTS: 
#'      1. linelist: the line list for observed individuals
#'      2. kinetics_pars: vector of named parameters for the viral kinetics model
#' OUTPUTS: 
#'      1. A tibble with the line list data and the viral load/ct/observed ct at the time of sampled
simulate_viral_loads_wrapper <- function(linelist,
                                         kinetics_pars){
  ## Control for changing standard deviation
  test_ages <- 0:1000
  t_switch <-  kinetics_pars["t_switch"] + kinetics_pars["desired_mode"] + kinetics_pars["tshift"]
  sd_mod <- rep(kinetics_pars["sd_mod"], max(test_ages))
  unmod_vec <- 1:min(t_switch,max(test_ages))
  sd_mod[unmod_vec] <- 1
  decrease_vec <- (t_switch+1):(t_switch+kinetics_pars["sd_mod_wane"])
  sd_mod[decrease_vec] <- 1 - ((1-kinetics_pars["sd_mod"])/kinetics_pars["sd_mod_wane"])*seq_len(kinetics_pars["sd_mod_wane"])
  vl_dat <- linelist %>% 
    ## Fix infection time for uninfected individuals
    mutate(infection_time = ifelse(is.na(infection_time),-100, infection_time)) %>%
    group_by(i) %>% 
    ## Pre-compute loss of detectability
    mutate(last_detectable_day = infection_time + ## From infection time
             rnbinom(n(), 1, prob=kinetics_pars["prob_detect"]) + ## How long until full clearance?
             kinetics_pars["tshift"] + kinetics_pars["desired_mode"] + kinetics_pars["t_switch"]) %>% ## With correct shift
    mutate(ct=virosolver::viral_load_func(kinetics_pars, sampled_time, FALSE, infection_time)) %>%
    ## If sampled after loss of detectability or not infected, then undetectable
    mutate(ct = ifelse(sampled_time > last_detectable_day, 1000, ct),
           ct = ifelse(is_infected == 0, 1000, ct),
           ct = ifelse(sampled_time < infection_time, 1000, ct)) %>%
    ## Convert to viral load value and add noise to Ct
    mutate(vl = ((kinetics_pars["intercept"]-ct)/log2(10)) + kinetics_pars["LOD"],
           days_since_infection = pmax(sampled_time - infection_time,-1),
           sd_used = ifelse(days_since_infection > 0, kinetics_pars["obs_sd"]*sd_mod[days_since_infection],kinetics_pars["obs_sd"]),
           ct_obs_sim = extraDistr::rgumbel(n(), ct, sd_used)) %>%
    ## Convert <LOD to LOD
    mutate(ct_obs = pmin(ct_obs_sim, kinetics_pars["intercept"])) %>%
    ungroup() %>%
    mutate(infection_time = ifelse(infection_time < 0, NA, infection_time))
  return(viral_loads=vl_dat)
}
