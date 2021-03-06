.datatable.aware <- TRUE
C_SEIR_model_rlsoda <- NULL
C_SEEIRR_model_rlsoda <- NULL
C_SEIR_switch_rlsoda <- NULL
C_initmodSEIR <- NULL
C_SEIR_model_lsoda <- NULL

#' @useDynLib virosolver
#' @importFrom Rcpp sourceCpp evalCpp
.onLoad <- function(...) {
  C_SEIR_model_rlsoda <<- getNativeSymbolInfo("SEIR_model_rlsoda", PACKAGE = "virosolver")
  C_SEEIRR_model_rlsoda <<- getNativeSymbolInfo("SEEIRR_model_rlsoda", PACKAGE = "virosolver")
  C_SEIR_switch_rlsoda <<- getNativeSymbolInfo("SEIR_switch_rlsoda", PACKAGE = "virosolver")
  
  C_SEIR_model_lsoda <<- getNativeSymbolInfo("SEIR_model_lsoda", PACKAGE = "virosolver")
  C_initmodSEIR <<- getNativeSymbolInfo("initmodSEIR", PACKAGE = "virosolver")
}

.onUnload <- function (libpath) {
  library.dynam.unload("virosolver", libpath)
}