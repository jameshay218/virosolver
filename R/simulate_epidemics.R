#' Wrapper for SEIR model simulation, deterministic or stochastic
#' 
#' Simulates a deterministic or stochastic SEIR model from model parameters
#' 
#'  @param N Population size. Defaults to 100000.
#'  @param times Vector of times to solve model over
#'  @param pars Named vector of SEIR model parameters
#'  @param version String either "ode" for the deterministic version or "odin" for stochastic
#'  @param switch_model If TRUE, uses the SEIR model with 2 switch points in transmission intensity
#'  @param beta_smooth spar parameter fed to smooth.spline, used to smooth the time-varying beta
#' 
#' @return Returns a list of 6 things: 
#'      1. Data frame of the SEIR solution
#'      2. Vector of absolute incidence per time point
#'      3. Overall probability of infection
#'      4. Plot of incidence over time
#'      5. Average growth rates over different time periods
#'      6. Plot of average growth rates over different time periods
#' @author James Hay, \email{jameshay218@@gmail.com}
#' @family simulation functions
#' 
#' @export
simulate_seir_wrapper <- function(N=100000, times, pars, version="ode",switch_model=FALSE,beta_smooth=0.8,
                                  convert_daily=TRUE){
  ## Choose which version of the model to solve
  if(version == "ode"){
    ## Deterministic SEIR model
    epidemic_process <- simulate_seir_process(pars,times,N=N,switch_model=switch_model)
    
    ## Widen solution and extract key variables
    res <- epidemic_process$seir_outputs %>% pivot_wider(names_from="variable",values_from="value")
    res <- res %>% rename(step=time)
    incidence <- epidemic_process$per_cap_incidence
    overall_prob <- epidemic_process$overall_prob_infection
  } else {
    ####################################################
    ## Stochastic model
    ####################################################
    if(switch_model){
      gamma1 <- 1/pars["infectious"]
      sigma1 <- 1/pars["incubation"]
      beta1 <- pars["R0_1"]*gamma1
      beta2 <- pars["R0_2"]*gamma1
      beta3 <- pars["R0_3"]*gamma1
      I0 <- ceiling(pars["I0"]*N)
      ## Odin stochastic SEIR model generator
      betas <- rep(beta3, length(times))
      betas[which(times < pars["t_switch2"])] <- beta2
      betas[which(times < pars["t_switch1"])] <- beta1
      betas <- smooth.spline(betas,spar=beta_smooth)$y
      seir <- seir_generator_interpolate$new(betat=times,betay=betas,sigma=sigma1,gamma=gamma1,S_ini=N-I0,I_ini=I0)
      
    } else {
      gamma1 <- 1/pars["infectious"]
      sigma1 <- 1/pars["incubation"]
      beta1 <- pars["R0"]*gamma1
      I0 <- ceiling(pars["I0"]*N)
      ## Odin stochastic SEIR model generator
      seir <- seir_generator$new(beta=beta1,sigma=sigma1,gamma=gamma1,S_ini=N-I0,I_ini=I0)
    }
    ## Solve model
    res <- seir$run(times)
    ## Make sure we get a simulation with an outbreak - keep trying until it takes off
    while(max(res[,"I"]) <= I0) res <- seir$run(times)
    res <- as.data.frame(res)
    
    ## Shift for start
    res$step <- res$step + floor(pars["t0"])
    
    ## Dummy rows from pre-seeding
    if(pars["t0"] > 0){
      if("incubation" %in% names(pars)){
        dummy_row <- data.frame("step"=0:(floor(unname(pars["t0"]))-1),"S"=N,"E"=0,"I"=0,"R"=0,"inc"=0)
      } else {
        dummy_row <- data.frame("step"=0:(floor(unname(pars["t0"]))-1),"S"=N,"I"=0,"R"=0,"inc"=0)
      }
      res <- bind_rows(dummy_row, res)
    }
    res <- res[res$step %in% times,]
    
    ## Get raw incidence and overall probability of infection
    incidence <- res$inc/N
    overall_prob <- max(res$R)/N
    
    if(!switch_model){
      res$beta <- pars["R0"]/pars["infectious"]
    }
    
    res$Rt <- (res$S/N) * res$beta * pars["infectious"]
  } 
  ## Get absolute incidence
  incidence <- incidence * N
  ## Reshape solution
  res_melted <- res %>% pivot_longer(-step)
  if(convert_daily){
    res_melted %>% mutate(step = floor(step)) %>% 
      group_by(name, step) %>% 
      summarize(value=mean(value)) %>%
      ungroup()
    
    sum_n_vector <- function (v, n) 
    {
      nv <- length(v)
      if (nv%%n) 
        v[ceiling(nv/n) * n] <- NA
      colSums(matrix(v, n), na.rm = TRUE)
    }
    
    incidence <- sum_n_vector(incidence, 1/diff(times)[1])*diff(times)[1]
    times <- unique(floor(times))
  }
  ## Compartment plot
  p_compartments <- res_melted %>% 
    dplyr::filter(name %in% c("S","E","I","R","cumulative_incidence")) %>%
    ggplot() + geom_line(aes(x=step,y=value,col=name)) +
    ylab("Per capita") +
    xlab("Date") +
    theme_bw() +
    theme(legend.position="top")
  
  ## Incidence plot
  p_inc <- ggplot(data.frame(x=times,y=incidence)) +
    geom_line(aes(x=x,y=y),col="red") +
    ylab("True incidence") +
    xlab("Date") +
    theme_bw()
  
  ## Rt plot
  p_rt <- res_melted %>% dplyr::filter(name == "Rt") %>%
    ggplot() +
    geom_line(aes(x=step,y=value),col="blue") +
    scale_y_continuous(limits=c(0,pars["R0"]+1)) +
    geom_hline(yintercept=1,linetype="dashed") +
    ylab("Rt") +
    xlab("Date") +
    theme_bw()
  
  ## Combine plots
  inc_plot <- p_compartments / p_inc / p_rt
  
  ## Get growth rates
  GR_daily <- log(incidence[2:length(incidence)]/incidence[1:(length(incidence)-1)])
  GR_daily <- ifelse(is.finite(GR_daily), GR_daily, NA)
  GR_daily_dat <- data.frame(t=times[2:length(times)],GR=GR_daily,ver="daily")
  
  ## Get average growth rate over different size windows
  lastdays <- seq(10,50,by=10)
  
  GR_all <- GR_daily_dat
  for(t in seq_along(lastdays)){
    GR_full <- NULL
    lastday <- lastdays[t]
    for (i in (lastday+1):length(times)) {
      end_index <- i-1
      start_index <- max(1, (i-lastday))
      GR_full <- c(GR_full,mean(GR_daily[start_index:end_index], na.rm=TRUE))
    }
    GR_full_dat <- data.frame(t=(lastday+1):length(times), GR=GR_full,ver=as.character(lastday))
    GR_all <- bind_rows(GR_all, GR_full_dat)
  }
  
  ## Get daily growth rate around the peak
  gr_crossover <- GR_all %>% dplyr::filter(ver == "daily") %>%
    mutate(abs_gr = abs(GR)) %>%
    dplyr::filter(abs_gr == min(abs_gr, na.rm=TRUE)) %>% pull(t)
  
  ## Average growth rates
  p_gr <- ggplot(GR_all) +
    geom_line(aes(x=t,y=GR,col=ver)) +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_vline(xintercept=gr_crossover,linetype="dotted")+
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    xlab("Date") +    
    ggtitle("Average growth rate over different windows") +
    theme_bw()
  
  ## Daily growth rate
  p_gr1 <- ggplot(GR_all %>% dplyr::filter(ver == "daily")) +
    geom_line(aes(x=t,y=GR),col="black") +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_vline(xintercept=gr_crossover,linetype="dotted")+
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    ggtitle("Daily growth rate") +
    xlab("Date") +    
    theme_bw()
  
  
  list(seir_outputs=res, 
       incidence=incidence,
       per_cap_incidence=incidence/N,
       overall_prob=overall_prob,
       plot=inc_plot, 
       growth_rates=GR_all,
       growth_rate_p =p_gr1/p_gr)
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
#' @export

simulate_seir_process <- function(pars, times, N=100000, switch_model=FALSE){
  if(switch_model){
    ## Pull parameters for SEIR model
    seir_pars <- c(pars["R0_1"]*(1/pars["infectious"]),pars["R0_2"]*(1/pars["infectious"]),pars["R0_3"]*(1/pars["infectious"]),
                   1/pars["incubation"],1/pars["infectious"],
                   pars["t_switch1"],pars["t_switch2"])
    
    ## Set up initial conditions.
    ## Note if population_n=1, then solves everything per capita
    init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0)
    
    ## Solve the SEIR model using the rlsoda package
    sol <- rlsoda::rlsoda(init, times, C_SEIR_switch_rlsoda, parms=seir_pars, dllname="virosolver",
                          deSolve_compatible = TRUE,return_time=TRUE,return_initial=TRUE,atol=1e-10,rtol=1e-10)
  } else {
    if("incubation" %in% names(pars)){
    ## Pull parameters for SEIR model
      seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
      ## Set up initial conditions.
      ## Note if population_n=1, then solves everything per capita
      # init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0)
      init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0)
      
        ## Solve the SEIR model using the rlsoda package
        #sol <- rlsoda::rlsoda(init, times, C_SEIR_model_rlsoda, parms=seir_pars, dllname="virosolver",
        #                      deSolve_compatible = TRUE,return_time=TRUE,return_initial=TRUE,atol=1e-10,rtol=1e-10)
        ## Solve the SEIR model using the lsoda package. lsoda runs about 4x slower than rlsoda, but
        ## the lsoda package is available on CRAN, making it more user-friendly.
        sol <- deSolve::ode(init, times, func="SEIR_model_lsoda",parms=seir_pars,
                            dllname="virosolver",initfunc="initmodSEIR",
                            nout=0, rtol=1e-10,atol=1e-10)
        dummy_row <- data.frame("time"=0:(floor(unname(pars["t0"]))-1),"S"=N,"E"=0,"I"=0,"R"=0,"cumu_exposed"=0,"Rt"=unname(pars["R0"])*N)
        use_colnames <- c("time","S","E","I","R","cumu_exposed")
        
    } else {
      sir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["infectious"])
      ## Set up initial conditions.
      ## Note if population_n=1, then solves everything per capita
      # init <- c((1-pars["I0"])*N,0,pars["I0"]*N,0,0)
      init <- c((1-pars["I0"])*N,pars["I0"]*N,0,0)
      
      ## Solve the SEIR model using the rlsoda package
      sol <- deSolve::ode(init, times, func="SIR_model_lsoda",parms=sir_pars,
                          dllname="virosolver",initfunc="initmodSIR",
                          nout=0, rtol=1e-10,atol=1e-10)
      dummy_row <- data.frame("time"=0:(floor(unname(pars["t0"]))-1),"S"=N,"I"=0,"R"=0,"cumu_exposed"=0,"Rt"=unname(pars["R0"])*N)
      use_colnames <- c("time","S","I","R","cumu_exposed")
    }
  }
  

  ## Convert to data frame and column names
  sol <- as.data.frame(sol)
  colnames(sol) <- use_colnames
  ## Get Rt
  sol$Rt <- (sol$S) * pars["R0"]
  
  ## Shift time for start
  sol$time <- sol$time + floor(pars["t0"])
  
  ## Dummy rows from pre-seeding
  if(pars["t0"] > 0){
    sol <- bind_rows(dummy_row, sol)
  }
  sol <- sol[sol$time %in% times,]
  
  ## Pull raw incidence in absolute numbers
  inc <- c(pars["I0"],diff(sol[,"cumu_exposed"]))
  inc <- pmax(0, inc)
  
  ## Get per capita incidence and prevalence
  per_cap_inc <- inc/N
  if("incubation" %in% names(pars)){ 
    per_cap_prev <- (sol$E + sol$I)/N
  } else {
    per_cap_prev <- (sol$I)/N
  }
  
  
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




## Adapted from http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(inc) <- n_SE
  
  ## Individual probabilities of transition:
  p_SE <- 1 - exp(-beta * I / N) # S to E
  p_EI <- 1 - exp(-sigma) # E to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_EI <- rbinom(E, p_EI)
  n_IR <- rbinom(I, p_IR)
  
  ## Total population size
  N <- S + E + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(inc) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  E_ini <- user(0)
  I_ini <- user(1)
  beta <- user(0.2)
  sigma <- user(0.15)
  gamma <- user(0.1)
  
}, verbose = FALSE)



## Adapted from http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator_switch <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(inc) <- n_SE
  
  ## Individual probabilities of transition:
  beta <- beta1*(step <= t_switch1) + beta2*(step > t_switch1 && step < t_switch2) + beta3*(step >= t_switch2)
  
  p_SE <- 1 - exp(-beta * I / N) # S to E
  p_EI <- 1 - exp(-sigma) # E to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_EI <- rbinom(E, p_EI)
  n_IR <- rbinom(I, p_IR)
  
  ## Total population size
  N <- S + E + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(inc) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  E_ini <- user(0)
  I_ini <- user(1)
  beta1 <- user(0.2)
  beta2 <- user(0.2)
  beta3 <- user(0.2)
  sigma <- user(0.15)
  gamma <- user(0.1)
  t_switch1 <- user(25)
  t_switch2 <- user(50)
  
  
}, verbose = FALSE)


## Adapted from http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator_interpolate <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(inc) <- n_SE
  
  ## Individual probabilities of transition:
  p_SE <- 1 - exp(-beta * I / N) # S to E
  p_EI <- 1 - exp(-sigma) # E to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_EI <- rbinom(E, p_EI)
  n_IR <- rbinom(I, p_IR)
  
  ## Total population size
  N <- S + E + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(inc) <- 0
  
  ## User defined parameters - default in parentheses:
  beta = interpolate(betat,betay,"constant")
  output(beta) = beta
  betat[] = user()
  betay[] = user()
  dim(betat) <- user()
  dim(betay) <- length(betat)
  
  S_ini <- user(1000)
  E_ini <- user(0)
  I_ini <- user(1)
  sigma <- user(0.15)
  gamma <- user(0.1)
  
  
}, verbose = FALSE)
