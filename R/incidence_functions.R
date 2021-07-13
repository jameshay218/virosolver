#' Wrapper function for the solveSEIRmodel_lsoda code, which is written in C. 
#' 
#' This function returns a vector of incidence values for given times 
#' according to SEIR model parameters. This function utilizes lsoda, 
#' a solver for ordinary differential equations. 
#' 
#' @param pars SEIR model parameters (dataframe)
#' @param times Time points for incidence calculation
#' 
#' @return Returns incidence values for time points passed in via the 'times' parameter.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' 
#' @export
solveSEIRModel_lsoda_wrapper <- function(pars, times){ 
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
  ## Initialize values for ODE system
  init <- c(1-pars["I0"],0,pars["I0"],0,0)
  ## Solve SEIR ordinary differential equations and calculate incidence for given time points
  inc <- c(rep(0, pars["t0"]),0,diff(
    deSolve::ode(init, times, func="SEIR_model_lsoda",parms=seir_pars, 
                 dllname="virosolver",initfunc="initmodSEIR",
                 nout=0, rtol=1e-6,atol=1e-6)[,6]))[1:length(times)]

  inc <- pmax(0, inc) ##FIX ME,  why restrict to positive values?  
  inc
}

#' Wrapper function for the solveSEIRmodel_rlsoda code, which is written in C. 
#' 
#' This function returns a vector of incidence values for given times 
#' according to SEIR model parameters. This function utilizes rlsoda, 
#' a solver for ordinary differential equations. 
#' 
#' @param pars SEIR model parameters (dataframe)
#' @param times Time points for incidence calculation
#' 
#' @return Returns incidence values for time points passed in via the 'times' parameter.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' 
#' @export
solveSEIRModel_rlsoda_wrapper <- function(pars, times){
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
  ## Initialize values for ODE system
  init <- c(1-pars["I0"],0,pars["I0"],0,0)
  ## Solve SEIR ordinary differential equations and calculate incidence for given time points
  inc <- c(rep(0, pars["t0"]),0,diff(rlsoda::rlsoda(init, times, C_SEIR_model_rlsoda, parms=seir_pars, dllname="virosolver",
                 deSolve_compatible = FALSE,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)[6,]))[1:length(times)]
  ## Removes negative values from  incidence vector
  inc <- pmax(0, inc)
  inc
}

#' Wrapper function for the solveSEIRswitch_rlsoda code, which is written in C. 
#' 
#' This function returns a vector of incidence values for given times 
#' according to SEIR model parameters. This function takes into account 
#' models in which R0 changes with time. 
#' 
#' This function utilizes rlsoda, a solver for ordinary differential equations. 
#' 
#' @param pars SEIR model parameters (dataframe)
#' @param times Time points for incidence calculation
#' 
#' @return Returns incidence values for time points passed in via the 'times' parameter.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' 
#' @export
solveSEIRswitch_rlsoda_wrapper <- function(pars, times){
  seir_pars <- c(pars["R0_1"]*(1/pars["infectious"]),
                 pars["R0_2"]*(1/pars["infectious"]),
                 pars["R0_3"]*(1/pars["infectious"]),
                 1/pars["incubation"],1/pars["infectious"],
                 pars["t_switch1"], pars["t_switch2"])
  ## Initialize values for ODE system
  init <- c(1-pars["I0"],0,pars["I0"],0,0)
  ## Solve SEIR ordinary differential equations and calculate incidence for given time points
  inc <- c(rep(0, pars["t0"]),0,diff(rlsoda::rlsoda(init, times, C_SEIR_switch_rlsoda, parms=seir_pars, dllname="virosolver",
                                                    deSolve_compatible = FALSE,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)[6,]))[1:length(times)]
  ## Removes negative values from  incidence vector
  inc <- pmax(0, inc)
  inc
}

#' Wrapper function for the solveSEIRRModel_rlsoda code, which is written in C. 
#' 
#' This function returns a vector of incidence values for given times 
#' according to SEIRR model parameters.
#' 
#' This function utilizes rlsoda, a solver for ordinary differential equations. 
#' 
#' @param pars SEIRR model parameters (dataframe)
#' @param times Time points for incidence calculation
#' 
#' @return Returns incidence values for time points passed in via the 'times' parameter.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' 
#' @export
solveSEEIRRModel_rlsoda_wrapper <- function(pars, times){
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["latent"],1/pars["incubation"],1/pars["infectious"],1/pars["recovery"])
  ## Initialize values for ODE system
  init <- c(1-pars["I0"],0,0,pars["I0"],0,0,0) ## FIXME : Why the additional 0's? seir_pars has 5, init has 7?
  ## Solve SEEIRR ordinary differential equations and calculate incidence for given time points
  inc <- c(rep(0, pars["t0"]),0,diff(rlsoda::rlsoda(init, times, C_SEEIRR_model_rlsoda, parms=seir_pars, dllname="virosolver",
                                                    deSolve_compatible = FALSE,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)[8,]))[1:length(times)]
  ## Removes negative values from  incidence vector
  inc <- pmax(0, inc)
  inc
}


#' This function returns a vector of detectable prevalence 
#' values for given times according to SEEIRR model parameters. 
#' 
#' This function utilizes C_SEEIRR_model_rlsoda, a function included 
#' in the virosolver that is written in C (code can be found in /src/)
#' and rlsoda, a solver for ordinary differential equations. 
#' 
#' @param pars SEEIRR model parameters (dataframe)
#' @param times Time points for prevalence calculation (FIXME: prevalence??)
#' 
#' @return Returns prevalence values for time points passed in via the 'times' parameter.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
detectable_SEEIRRModel <- function(pars, times){
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["latent"],1/pars["incubation"],1/pars["infectious"],1/pars["recovery"])
  ## Initialize values for ODE system
  init <- c(1-pars["I0"],0,0,pars["I0"],0,0,0)
  ## Solve SEEIRR ordinary differential equations and calculate incidence for given time points
  y <- rlsoda::rlsoda(init, times, C_SEEIRR_model_rlsoda, parms=seir_pars, dllname="virosolver",
                      deSolve_compatible = FALSE,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)
  detectable_prev <- colSums(y[c(4,5,6),])
  detectable_prev <- c(rep(0, pars["t0"]),detectable_prev)
  ## Removes negative values from prevalence vectors
  pmax(0, detectable_prev)
}


#' This function returns a dataframe containing the values solved for  
#' using the SEEIRR differential equation set
#' 
#' @param pars SEEIRR model parameters (dataframe)
#' @param times Time points for prevalence calculation 
#' 
#' @return Returns dataframe of variables solved for using SEEIRR ODEs
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
SEEIRRModel <- function(pars, times){
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["latent"],1/pars["incubation"],1/pars["infectious"],1/pars["recovery"])
  init <- c(1-pars["I0"],0,0,pars["I0"],0,0,0)
  ## Solves for SEEIRR ODEs and converts output into dataframe
  y <- rlsoda::rlsoda(init, times, C_SEEIRR_model_rlsoda, parms=seir_pars, dllname="virosolver",
                      deSolve_compatible = TRUE,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)
  y <- as.data.frame(y)
  ## Shifts times in relation to initial time 
  y$time <- y$time + floor(pars["t0"])
  addition <- as.data.frame(matrix(0, ncol=ncol(y),nrow=floor(pars["t0"])))
  colnames(addition) <- colnames(y)
  addition$time <- 0:(floor(pars["t0"])-1)
  ## Attaches time column to dataframe
  y_end <- rbind(addition,y) ## FIXME: running into error on this step --- come back and double-check. 
  colnames(y_end) <- colnames(y)
  y_end
}

#' This function returns incidence for given times 
#' using the exponential growth model. 
#' 
#' @param pars Exponential growth model parameters
#' @param times Time points for incidence calculation 
#' 
#' @return Returns incidence for given time points
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
exponential_growth_model <- function(pars, times){
  overall_prob <- pars["overall_prob"]
  beta <- pars["beta"]

  y <- exp(beta*times)
  scale <- sum(y)/overall_prob
  prob_infection_tmp <- y/scale
  prob_infection_tmp
}

#' This function returns incidence for given times 
#' using the gaussian process model. 
#' 
#' @param pars Gaussan process model parameters
#' @param times Time points for incidence calculation 
#' 
#' @return Returns incidence for given time points
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
gaussian_process_model <- function(pars, times){
  par_names <- names(pars)
  use_names <- c("prob",paste0("prob.",1:length(times)))
  overall_prob <- pars["overall_prob"]
  k <- pars[which(par_names%in%use_names)]

  mat <- matrix(rep(times, each=length(times)),ncol=length(times))
  t_dist <- abs(apply(mat, 2, function(x) x-times))
  mus <- rep(0, length(times))
  nu <- pars["nu"]
  rho <- pars["rho"]
  K <- nu^2 * exp(-rho^2 * t_dist^2)
  diag(K) <- diag(K) + 0.01
  L_K <- t(chol(K))
  k1 <- (L_K %*% k)[,1]
  ps <- 1/(1+exp(-k1))
  ps <- ps/sum(ps)
  ## Calculates incidence for given times 
  prob_infection_tmp <- ps*overall_prob
  prob_infection_tmp
}

#' Given a vector of daily infection probabilities, 
#' this function converts this to the input expected
#' by the Gaussian process model
#' 
#' @param desired_probs FIX ME 
#' @param pars Exponential growth model parameters
#' @param times Time points for incidence calculation 
#' 
#' @return Returns daily infection probabilities in the format 
#' expected by the Gaussian process model
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
reverse_gp_model <- function(desired_probs, pars, times){ ## FIXME: review upderpinnings of this function
  mat <- matrix(rep(times, each=length(times)),ncol=length(times))
  t_dist <- abs(apply(mat, 2, function(x) x-times))
  nu <- pars["nu"]
  rho <- pars["rho"]
  K <- nu^2 * exp(-rho^2 * t_dist^2)
  diag(K) <- diag(K) + 0.01
  L_K <- t(chol(K)) ## FIXME: where is "t" initialized?


  ps <- sum(desired_probs)*desired_probs + 0.000000001
  ps <- pars["overall_prob"]*ps/sum(ps)
  k1 <- -log((1/ps) -1)

  k_solve <- (k1 %*% solve(t(L_K)))[1,]
  k_solve
}


#' This function solves the SEIR ODEs and returns
#' a matrix of class deSolve. This function utilizes
#' ODE solver deSolve, which intrinsically uses lsoda. 
#' 
#' @param init Initial conditions for system of equations 
#' @param pars SEIR model parameters
#' @param times Time points for incidence calculation 
#' 
#' @return Returns matrix of class deSolve containing solutions for SEIR ODEs 
#' for each time in ts 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
solveSEIRModel_lsoda <- function(ts, init, pars,compatible=FALSE){
  pars <- pars[c("beta","sigma","gamma")] ## FIXME: which parameters are represented by beta, sigma and gamma? 
  deSolve::ode(init, ts, func="SEIR_model_lsoda",parms=pars,
               dllname="virosolver",initfunc="initmodSEIR",
               nout=0, rtol=1e-6,atol=1e-6)
}

#' This function solves the SEIR ODEs for each time in ts and returns 
#' a matrix of values solved for the time points specified in 'ts'. 
#' 
#' @param ts Time points for incidence calculation 
#' @param init Initial conditions for system of equations 
#' @param pars SEIR model parameters
#' @param compatible A boolean variable indicating if
#' results should be transposed in deSolve compatible format.
#' Set to FALSE by default. 
#' 
#' @return Returns matrix of class deSolve containing solutions for SEIR ODEs 
#' for each time in ts 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
solveSEIRModel_rlsoda <- function(ts, init, pars,compatible=FALSE){
  pars <- pars[c("beta","sigma","gamma")]
  ## Solves for SEIR ODEs
  rlsoda::rlsoda(init, ts, C_SEIR_model_rlsoda, parms=pars, dllname="virosolver",
                 deSolve_compatible = compatible,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)
}

#' This function solves the SEIR ODEs for each time in ts and returns 
#' a matrix of values solved for the time points specified in 'ts'. 
#' This function takes into account models in which R0 changes with time.
#' 
#' @param ts Time points for incidence calculation 
#' @param init Initial conditions for system of equations 
#' @param pars SEIR model parameters
#' @param compatible A boolean variable indicating if
#' results should be transposed in deSolve compatible format.
#' Set to FALSE by default. 
#' 
#' @return Returns matrix of class deSolve containing solutions for SEIR ODEs 
#' for each time in ts 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
solveSEIRswitch_rlsoda <- function(ts, init, pars,compatible=FALSE){
  pars <- pars[c("beta1","beta2","beta3","sigma","gamma","t_switch1","t_switch2")]
  print(pars)
  ## Solves for SEIR ODEs
  rlsoda::rlsoda(init, ts, C_SEIR_switch_rlsoda, parms=pars, dllname="virosolver",
                 deSolve_compatible = compatible,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)
}

#' This function solves the SEEIRR ODEs for each time in ts and returns 
#' a matrix of values solved for the time points specified in 'ts'. 
#' 
#' @param ts Time points for incidence calculation 
#' @param init Initial conditions for system of equations 
#' @param pars SEEIRR model parameters
#' @param compatible A boolean variable indicating if
#' results should be transposed in deSolve compatible format.
#' Set to FALSE by default. 
#' 
#' @return Returns matrix of class deSolve containing solutions for SEEIRR ODEs 
#' for each time in ts 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
solveSEEIRRModel_rlsoda <- function(ts, init, pars,compatible=FALSE){
  pars <- pars[c("beta","sigma","alpha","gamma","omega")]
  ## Solves for SEIR ODEs
  rlsoda::rlsoda(init, ts, C_SEEIRR_model_rlsoda, parms=pars, dllname="virosolver",
                 deSolve_compatible = compatible,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)
}

#' This function solves the SEIR ODEs for each time in ts and returns 
#' a vector of detectable prevalence values for each time point given in
#' the 'times' parameter
#' 
#' @param pars SEIR model parameters
#' @param times Time points for incidence calculation 
#' 
#' @return 
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family incidence functions
#' 
#' @examples FIX ME
#' @export
detectable_SEIRModel <- function(pars, times){
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
  init <- c(1-pars["I0"],0,pars["I0"],0)

  ## Solves for SEIR ODEs
  y <- rlsoda::rlsoda(init, times, C_SEIR_model_rlsoda, parms=seir_pars, dllname="virosolver",
                 deSolve_compatible = FALSE,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)

  detectable_prev <- y[4,]
  detectable_prev <- c(rep(0, pars["t0"]),detectable_prev)
  pmax(0, detectable_prev)
}
