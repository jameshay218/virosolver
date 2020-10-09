#' @export
solveSEIRModel_rlsoda_wrapper <- function(pars, times){
  seir_pars <- c(pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
  init <- c(1-pars["I0"],0,pars["I0"],0,0)
  c(rep(0, pars["t0"]),0,diff(rlsoda::rlsoda(init, times, C_SEIR_model_rlsoda, parms=seir_pars, dllname="virosolver",
                 deSolve_compatible = FALSE,return_time=TRUE,return_initial=TRUE,atol=1e-5,rtol=1e-5)[6,]))[1:length(times)]
}


#' @export
solveSEIRModel_rlsoda <- function(ts, init, pars,compatible=FALSE){
  pars <- pars[c("beta","sigma","gamma")]
  rlsoda::rlsoda(init, ts, C_SEIR_model_rlsoda, parms=pars, dllname="virosolver",
                 deSolve_compatible = compatible,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)
}
