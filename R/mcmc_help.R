#' Function to generate random starting parameter values that 
#' return a finite likelihood.
#' 
#' This function is used to generate random starting values 
#' for the MCMC procedure.
#' 
#' @param partab Dataframe containing parameters corresponding
#' to underlying incidence model (e.g. SEIR, SEEIRR).
#' @param obs_dat Dataframe containing observed Ct values.
#' @param CREATE_POSTERIOR_FUNC Creates the posterior function 
#' used in the MCMC framework for detectable Ct values.
#' @param INCIDENCE_FUNC Function that expects a vector of named parameters 
#' and returns a vector of daily incidence.
#' @param PRIOR_FUNC Function that returns the log prior probability for a 
#' given vector of parameter values given the prior means and sds.
#' @param use_pos Boolean variable indicating if only positive Ct values
#' i.e, those below a specified threshold,  should be used. Set to FALSE 
#' by default.
#' 
#' @return Dataframe containing random starting parameter values that return 
#' a finite likelihood.
#' 
#' @author James Hay, \email{jhay@@hsph.harvard.edu}
#' @family mcmc help
#' 
#' @examples FIX ME
#' 
#' @export
generate_viable_start_pars <- function(parTab,
                                       obs_dat,
                                       CREATE_POSTERIOR_FUNC,
                                       INCIDENCE_FUNC,
                                       PRIOR_FUNC,
                                       use_pos=FALSE,
                                       ...){
  f <- CREATE_POSTERIOR_FUNC(parTab, obs_dat,
                             INCIDENCE_FUNC=INCIDENCE_FUNC,
                             PRIOR_FUNC=PRIOR_FUNC,
                             use_pos=use_pos,
                             ...)
  ## Generates random values between lower and upper start
  ##   defined in parTab
  startTab <- lazymcmc::generate_start_tab(parTab)
  lik <- f(startTab$values)
  ## Ensures that startTab parameters return finite likelihood
  while(!is.finite(lik)){
    startTab <- lazymcmc::generate_start_tab(startTab)
    lik <- f(startTab$values)
  }
  startTab
}
