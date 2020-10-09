#' @export
create_posterior_func <- function(parTab,
                                  data,
                                  PRIOR_FUNC=NULL,
                                  INCIDENCE_FUNC=NULL,
                                  solve_ver="likelihood",
                                  solve_likelihood=TRUE) {
  par_names <- parTab$names
  times <- 0:max(data$t)
  ages <- 1:max(data$t)
  obs_times <- unique(data$t)

  f <- function(pars){
    names(pars) <- par_names
    prob_infection_tmp <- INCIDENCE_FUNC(pars, times)

    if(solve_ver == "likelihood"){
      lik <- 0
      if(solve_likelihood){
        lik <- sum(likelihood_cpp_wrapper(data, ages, pars, prob_infection_tmp))
      }

      if(!is.null(PRIOR_FUNC)){
        prior <- PRIOR_FUNC(pars)
        lik <- lik+prior
      }
      return(lik)
    } else {
      preds <- pred_dist_wrapper(seq(0,40,by=1),obs_times,ages,pars,prob_infection_tmp)
      return(preds)
    }
  }
  f
}
