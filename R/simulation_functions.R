#' Simulate viral loads
#' 
#' Takes infection times, the times over which to solve the model, and 
#' the viral kinetics parameters and returns a tibble containing true and observed viral loads.
#' 
#' @param infection_times Vector containing infection times. 
#' @param solve_times Vector of times at which individuals can be reported.
#' @param kinetics_pars Vector of named parameters for the viral kinetics model.
#' @param vl_func_use Determines which viral load function to use. 
#' The default is viral_load_func.
#' @param convert_vl Convert to viral load value. FALSE by default.
#' @param add_noise Add noise to Ct values? NULL by default.
#' 
#' @return Returns a tibble containing true and observed viral loads.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family simulation functions
#' 
#' @examples FIX ME
#' 
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
  ## How many days do you remain detectable after hitting the switch point?
  days_still_detectable <- rnbinom(n, size=1,prob=kinetics_pars["prob_detect"])

  for(i in seq_along(infection_times)){
    ## %% is modular division used here as a progress indicator. Generates a message with
    ## i every 1000th iteration
    if(i %% 1000 == 0) message(i)
    if(infection_times[i] > 0){
      pars <- kinetics_pars
      mod_probs <- rep(1, length(solve_times))
      ## Calling the viral load function to get Ct values
      ct <- vl_func_use(pars, solve_times, FALSE, infection_times[i])
      ## Additional way that individuals can become undetectable
      ## How long to wait until undetectable?
      mod_probs[which(solve_times >= kinetics_pars["tshift"] +
                        kinetics_pars["desired_mode"] +
                        kinetics_pars["t_switch"] +
                        infection_times[i] +
                        days_still_detectable[i])] <- 0
      ct[ct > kinetics_pars["true_0"]] <- kinetics_pars["true_0"]
      ct[which(solve_times < infection_times[i])] <- 1000
      ct <- ct * mod_probs
      ct[which(mod_probs == 0)] <- 1000
      viral_loads[i,] <- ct
    }
  }
  colnames(viral_loads) <- solve_times

  ## Convert to viral load value
  if(convert_vl) {
    viral_loads <- ((kinetics_pars["intercept"]-viral_loads)/log2(10)) + kinetics_pars["LOD"]
    true_viral_loads <- viral_loads
    ## If some of the viral loads are less than the limit of detection, re-assign these values
    ## to the limit of detection
    true_viral_loads[true_viral_loads < kinetics_pars["LOD"]] <- kinetics_pars["LOD"]
  } else {
    true_viral_loads <- viral_loads
    ## If some of the viral loads are greater than the intercept parameter, re-assign these
    ## values to equal the intercept (the maximum number of cycles run on the PCR machine, always
    ## assumed to be 40 in the simulations).
    true_viral_loads[true_viral_loads > kinetics_pars["intercept"]] <- kinetics_pars["intercept"]
  }
  
  ## Add noise to Ct values
  if(!is.null(add_noise)){
    observed_viral_loads <- t(apply(viral_loads, 1, function(vl) add_noise(length(vl),vl, kinetics_pars["obs_sd"])))
    colnames(observed_viral_loads) <- solve_times
    if(convert_vl) {
      ## If some of the viral loads are less than the limit of detection, re-assign these values
      ## to the limit of detection
      observed_viral_loads[observed_viral_loads < kinetics_pars["LOD"]] <- kinetics_pars["LOD"]
    } else {
      ## If some of the viral loads are greater than the intercept parameter, re-assign these
      ## values to equal the intercept (the maximum number of cycles run on the PCR machine)
      observed_viral_loads[observed_viral_loads > kinetics_pars["intercept"]] <- kinetics_pars["intercept"]
    }
  } else {
    observed_viral_loads <- true_viral_loads
  }
  
  ## Transform wide-format data and into long-format data and add column names
  true_viral_loads_melted <- reshape2::melt(true_viral_loads)
  colnames(true_viral_loads_melted) <- c("i","t","true")
  obs_viral_loads_melted <- reshape2::melt(observed_viral_loads)
  colnames(obs_viral_loads_melted) <- c("i","t","obs")

  # Adds columns from y (observed viral loads) to x (true viral loads), including all rows in x 
  combined_dat <- left_join(true_viral_loads_melted,obs_viral_loads_melted) %>% as_tibble

  return(combined_dat)
}

#' Simulate infection times
#' 
#' Infection times are simulated using the probability of infection and sample size n.
#' 
#' @param n Sample size. Must be an integer.
#' @param prob_infection A vector containing probabilities of infection.
#' @param overall_prob The overall probability of infection. NULL by default.
#' 
#' @return Returns a vector of infection times.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family simulation functions
#' 
#' @examples FIX ME
#' 
#'@export

simulate_infection_times <- function(n, prob_infection, overall_prob=NULL){
  
  ## Sum the probability of infection to get overall probability of infection
  if(is.null(overall_prob)){
    overall_prob <- sum(prob_infection)
  }
  
  ## Scale the probability infection by the overall probability of infection
  scaled_prob<- prob_infection/sum(prob_infection)
  
  are_infected <- numeric(n)
  infection_times <- numeric(n)
  
  ## For each sample n, simulate infection from a binomial distribution using the 
  ## overall probability of infection. 
  for(i in 1:n){
    infection <- rbinom(1,1, overall_prob)
    are_infected[i] <- infection
    ## If infected, sample from the probabilities of infection using the scaled
    ## probability of infection
    if(infection == 1){
      t_inf <- sample(1:length(prob_infection), 1, prob=scaled_prob)
      infection_times[i] <- t_inf
    } else {
      ## If not infected, assign -1 to infection times
      infection_times[i] <- -1
    }
  }
  return(infection_times)
}

#' Simulate full line list data
#'
#' Takes a vector of infection incidence (absolute numbers) and returns a tibble with line list data for all
#' individuals in the population. Infected individuals are assigned symptom onset, incubation 
#' periods (if applicable) and confirmation delays from log normal and gamma distributions respectively.
#'
#' @param incidence A vector of infection incidence (absolute numbers).
#' @param times Calendar time (i.e. days since the start of the epidemic).
#' @param symp_frac Fraction of the population that is symptomatic. Defaults to 0.35.
#' @param population_n Size of the population. Defaults to 100,000. Must match the N used
#' in the simulate_seir_process function.
#' @param incu_period_par1 The first parameter (meanLog) associated with the log-normal incubation period.
#' Defaults to 1.621.
#' @param incu_period_par2 The second parameter (sdLog) associated with the log-normal incubation period. 
#' Defaults to 0.418.
#' @param conf_delay_par1 The first parameter (shape) associated with the discretized gamma confirmation delay. 
#' Defaults to 5.
#' @param conf_delay_par2 The second parameter (rate) associated with the discretized gamma confirmation delay.
#' Defaults to 2. 
#'
#' @return A tibble with line list data for all individuals in the population.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family simulation functions
#' 
#' @examples FIX ME
#'
#' @export

simulate_observations_wrapper <- function(
  incidence, times, symp_frac=0.35,
  population_n=100000,
  incu_period_par1=1.621,incu_period_par2=0.418,
  conf_delay_par1=5,conf_delay_par2=2){
  
  ## Number of individuals who are not infected in the population
  not_infected <- population_n-sum(incidence)
  ## Create a tibble with infection times and incidence
  inc_dat <- tibble(infection_time=c(NA,times),inc=c(not_infected,incidence))
  ## Duplicates rows in inc_dat according to inc
  inc_dat <- inc_dat %>% uncount(inc)

  inc_dat <- inc_dat %>%
    mutate(i=1:n()) %>%
    ## If infection time is missing, assign 0 to is_infected, otherwise assign 1
    mutate(is_infected = ifelse(is.na(infection_time), 0, 1)) %>%
    ## Is this individual going to be symptomatic? Uses a binomial distribution
    mutate(is_symp=ifelse(is_infected, rbinom(n(), 1, symp_frac), 0)) %>%
    ## Symptom onset time from a log normal distribution
    mutate(incu_period=ifelse(is_infected & is_symp, rlnorm(n(), incu_period_par1, incu_period_par2), NA),
           onset_time=infection_time+floor(incu_period)) %>%
    ## Confirmation time from a gamma distribution
    mutate(confirmation_delay=extraDistr::rdgamma(n(),shape=conf_delay_par1,rate=conf_delay_par2))
  inc_dat
}

#' Simulate SEIR model
#' 
#' Simulates a deterministic SEIR model from model parameters, times to solve over, and 
#' population size.
#' 
#' @param pars SEIR model parameters.
#' @param times Times over which the model is solved.
#' @param N Population size. Defaults to 100000.
#' 
#' @return Returns a list of 7 things: 
#' 1. Plot of all SEIR compartments over time
#' 2. Plot of incidence and prevalence over time
#' 3. Solution of ordinary differential equation
#' 4. Absolute incidence per time point (raw incidence)
#' 5. Per capita incidence per time point
#' 6. Per capita prevalence (compartments E+I) per time point
#' 7. Overall probability of infection
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family simulation functions
#' 
#' @examples FIX ME
#' 
#' @export

simulate_seir_process <- function(pars, times, N=100000){
  ## Pull parameters for SEIR model
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
  ## Set up initial conditions.
  ## Note if population_n=1, then solves everything per capita
  init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0,0)
  
  ## Solve the SEIR model using the rlsoda package
  #sol <- rlsoda::rlsoda(init, times, C_SEIR_model_rlsoda, parms=seir_pars, dllname="virosolver",
  #                      deSolve_compatible = TRUE,return_time=TRUE,return_initial=TRUE,atol=1e-10,rtol=1e-10)
  
  ## Solve the SEIR model using the lsoda package. lsoda runs about 4x slower than rlsoda, but
  ## the lsoda package is available on CRAN, making it more user-friendly.
  sol <- deSolve::ode(init, times, func="SEIR_model_lsoda",parms=seir_pars,
                      dllname="virosolver",initfunc="initmodSEIR",
                      nout=0, rtol=1e-6,atol=1e-6)
  
  ## Convert to data frame and column names
  sol <- as.data.frame(sol)
  colnames(sol) <- c("time","S","E","I","R","cumu_exposed","cumu_infected")
  ## Get Rt
  sol$Rt <- (sol$S) * pars["R0"]
  
  ## Shift time for start
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
  
  ## Melt solution (wide to long format) and get per capita
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
    ylab("Per capita incidence (red)\n and prevalence (blue)") +
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

#' Simulate reporting
#' 
#' Line list data are subset by testing strategy. There are four options:
#' 1. Sample a random fraction of the population if the only argument is frac_report.
#' 2. Sample some random fraction of the population at a subset of time points, specified 
#' by timevarying_prob.
#' 3. Observe symptomatic individuals with some fixed probability, frac_report if symptomatic 
#' is TRUE.
#' 4. Observe symptomatic individuals with some time-varying probability, timevarying_prob, 
#' if symptomatic is TRUE.
#' 
#' @param individuals The full line list from the simulation, returned by 
#' virosolver::simulate_observations_wrapper.
#' @param solve_times Vector of times at which individuals can be reported.
#' @param frac_report The overall fraction/probability of individuals who are reported. 
#' Defaults to 1.
#' @param timevarying_prob A tibble with variables t and prob. This gives the probability of 
#' being reported on day t. NULL by default.
#' @param symptomatic If TRUE, then individuals are reported after developing symptoms. 
#' If FALSE, then we take a random cross-section. Defaults to FALSE.
#' 
#' @return Returns a list of 3 things: 
#' 1. A tibble with line list data for individuals who were observed.
#' 2. A plot of incidence for both observed individuals and the entire simulated population.
#' 3. Plot growth rate of cases/infections in the entire population and observed population.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family simulation functions
#' 
#' @examples FIX ME
#' 
#' @export

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
    ## Base case, individuals are sampled at random
    if(!symptomatic){
      ## Sample a fraction of the overall line list and assign each individual a sample time (at random)
      ## and a confirmation time (with confirmation delay)
      sampled_individuals <- sampled_individuals %>% 
        sample_frac(frac_report) %>%
        group_by(i) %>%
        mutate(sampled_time=sample(solve_times,n()),
               ## We sample, but then have to wait for a result
               confirmed_time=sampled_time+confirmation_delay) %>%
        ungroup()
    } else {
      ## Symptomatic based surveillance. Observe individuals based on symptom onset date
      ## Subset by symptomatic individuals. Observe some fraction of these at some delay post symptom onset
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
      
    } else {
      ## This is quite different - if you have symptom onset on day t, there is a probability that you will be observed
      ## by symptomatic based surveillance. Observe individuals based on symptom onset date
      sampled_individuals <- sampled_individuals %>% 
        filter(is_symp==1) %>% ## Subset for symptomatic
        left_join(timevarying_prob %>% rename(onset_time=t) %>% ## Join with time-varying reporting rate table
                    dplyr::select(-ver), by="onset_time") %>%
        group_by(i) 
      
      ## Quicker to vectorize
      ## Assign individuals as detected or not
      sampled_individuals$is_reported <- rbinom(nrow(sampled_individuals), 1, sampled_individuals$prob) 
      
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
  
  ## Group entire dataset
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
  
  ## Plot growth rate of cases/infections in the entire population and observed population
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

#' Simulate observed Ct values for the line list dataset 
#' 
#' This differs from virosolver::simulate_viral_loads,
#' as this function only solves the viral kinetics model for the observation time.
#' 
#' @param linelist The line list for observed individuals.
#' @param kinetics_pars Vector of named parameters for the viral kinetics model.
#' 
#' @return A tibble with the line list data and the viral load/ct/observed ct at the time of 
#' sample collection.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family simulation functions
#' 
#' @examples FIX ME
#' 
#' @export

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
           days_since_infection = pmax(floor(incu_period + confirmation_delay),-1),
           sd_used = ifelse(days_since_infection > 0, kinetics_pars["obs_sd"]*sd_mod[days_since_infection],kinetics_pars["obs_sd"]),
           ct_obs_sim = round(extraDistr::rgumbel(n(), ct, sd_used),2)) %>%
    ## If some of the viral loads are greater than the intercept parameter, re-assign these
    ## values to equal the intercept (the maximum number of cycles run on the PCR machine, always
    ## assumed to be 40 in the simulations).
    mutate(ct_obs = pmin(ct_obs_sim, kinetics_pars["intercept"])) %>%
    ungroup() %>%
    mutate(infection_time = ifelse(infection_time < 0, NA, infection_time))
  
  return(viral_loads=vl_dat)
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


#' Simulate N Ct values under symptom-based surveillance
#' and returns a dataframe of Ct values and times since infection. 
#'  
#' @param ages Vector of ages (days since infection).
#' @param kinetics_pars Vector of named parameters for the viral kinetics model.
#' @param incu_par1 logMean of the log-normal incubation period distribution
#' @param incu_par2 logSD of the log-normal incubation period distribution
#' @param sampling_par1 shape of the discretized gamma sampling delay distribution 
#' @param sampling_par2 rate of the discretized gamma sampling delay distribution 
#' @param N Number of observations to randomly generate overall;
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
simulate_viral_loads_example_symptoms <- function(ages, kinetics_pars,
                                                  incu_par1=1.621,incu_par2=0.418,
                                                  sampling_par1=5,sampling_par2=1,
                                                  N=100){
  t_switch <-  kinetics_pars["t_switch"] + kinetics_pars["desired_mode"] + kinetics_pars["tshift"]
  sd_mod <- rep(kinetics_pars["sd_mod"], max(ages))
  unmod_vec <- 1:min(t_switch,max(ages))
  sd_mod[unmod_vec] <- 1
  decrease_vec <- (t_switch+1):(t_switch+kinetics_pars["sd_mod_wane"])
  ## For the next sd_mod_wane days, variance about Ct trajectory decreases linearly
  sd_mod[decrease_vec] <- 1 - ((1-kinetics_pars["sd_mod"])/kinetics_pars["sd_mod_wane"])*seq_len(kinetics_pars["sd_mod_wane"])
  
  modal_cts <- virosolver::viral_load_func(kinetics_pars, ages, FALSE, 0)
  
  ## Simulate an incubation period, sampling delay and time of undetectable for each individual
  ## For each value we're going to simulate, pre-determine if it will still be detectable or not
  detectable_statuses <- floor(rnbinom(N, 1, prob=kinetics_pars["prob_detect"]) + 
    kinetics_pars["tshift"] + kinetics_pars["desired_mode"] + kinetics_pars["t_switch"])
  incu_periods <- floor(rlnorm(N, incu_par1,incu_par2))
  sampling_delays <- extraDistr::rdgamma(N, sampling_par1, sampling_par2)
  days_since_infection <- incu_periods + sampling_delays
  
  sim_dat <- tibble(i=1:N, day_undetectable=detectable_statuses, incu_period=incu_periods,sampling_delay=sampling_delays,days_since_infection=days_since_infection)
  
  sim_dat <- sim_dat %>% group_by(i) %>%
    mutate(modal_ct=modal_cts[days_since_infection+1]) %>%
    mutate(ct_obs=extraDistr::rgumbel(n(), modal_ct, kinetics_pars["obs_sd"]*sd_mod[days_since_infection+1])) %>%
    mutate(ct_obs=pmin(ct_obs, kinetics_pars["intercept"])) %>%
    rename(age=days_since_infection)
  
  return(sim_dat)
}

