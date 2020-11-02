#' @export
create_posterior_func <- function(parTab,
                                  data,
                                  PRIOR_FUNC=NULL,
                                  INCIDENCE_FUNC=NULL,
                                  solve_ver="likelihood",
                                  solve_likelihood=TRUE,
                                  use_pos=FALSE,
                                  ...) {
  par_names <- parTab$names
  times <- 0:max(data$t)
  ages <- 1:max(data$t)
  obs_times <- unique(data$t)

  f <- function(pars){
    names(pars) <- par_names
    prob_infection_tmp <- INCIDENCE_FUNC(pars, times)
    #prob_infection_tmp[prob_infection_tmp < pars["I0"]] <- pars["I0"]
    if(solve_ver == "likelihood"){
      lik <- 0
      if(solve_likelihood){
        lik <- sum(likelihood_cpp_wrapper(data, ages, pars, prob_infection_tmp,use_pos))
      }

      if(!is.null(PRIOR_FUNC)){
        prior <- PRIOR_FUNC(pars, ...)
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

#' @export
create_post_func_seeirr <- function(parTab, data, ts, INCIDENCE_FUNC=detectable_SEEIRRModel, PRIOR_FUNC=NULL,ver="likelihood"){
  par_names <- parTab$names
  test_times <- unique(data$date)
  observed_prev <- data$POS
  N_obs <- data$POS + data$NEG

  f <- function(pars){
    names(pars) <- par_names
    prev <- INCIDENCE_FUNC(pars, ts)
    if(ver == "likelihood"){
      prev <- prev[which(ts %in% test_times)]
      lik <- sum(dbinom(observed_prev, N_obs, prev,log=TRUE))
      if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
      return(lik)
    } else {
      return(prev)
    }
  }
  f
}
#' @export
create_post_func_seeirr_combined <- function(parTab, data, ts, INCIDENCE_FUNC=detectable_SEEIRRModel,
                                             PRIOR_FUNC=NULL,ver="likelihood"){
  par_names <- parTab$names
  unique_locs <- unique(data$location)
  f <- function(pars){
    names(pars) <- par_names
    all_lik <- numeric(length(unique_locs))
    all_prev <- NULL

    for(i in seq_along(unique_locs)){
      loc <- unique_locs[i]
      tmp_dat <- data %>% filter(location == loc)
      test_times <- unique(tmp_dat$date)
      observed_prev <- tmp_dat$POS
      N_obs <- tmp_dat$POS + tmp_dat$NEG

      tmp_pars <- c(pars[paste0("R0_",i)], pars[paste0("t0_",i)],
                    pars["infectious"],pars["latent"],pars["incubation"],pars["recovery"],pars["I0"])
      names(tmp_pars) <- c("R0","t0","infectious","latent","incubation","recovery","I0")

      prev <- INCIDENCE_FUNC(tmp_pars, ts)[1:length(ts)]
      prev_subset <- prev[which(ts %in% test_times)]
      lik <- sum(dbinom(observed_prev, N_obs, prev_subset,log=TRUE))
      if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
      all_lik[i] <- lik
      if(ver == "model"){
        all_prev[[i]] <- tibble(prev=prev,loc=loc,t=ts)
      }
    }
    if(ver == "likelihood"){
      return(sum(all_lik))
    }
    if(ver == "model") {
      all_prev <- do.call("bind_rows", all_prev)
      return(all_prev)
    }
    return(all_lik)
  }
  f
}
