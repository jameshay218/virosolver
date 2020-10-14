.datatable.aware <- TRUE
C_SEIR_model_rlsoda <- NULL
C_SEEIRR_model_rlsoda <- NULL
#' @useDynLib virosolver
#' @importFrom Rcpp sourceCpp
.onLoad <- function(...) {
  C_SEIR_model_rlsoda <<- getNativeSymbolInfo("SEIR_model_rlsoda", PACKAGE = "virosolver")
  C_SEEIRR_model_rlsoda <<- getNativeSymbolInfo("SEEIRR_model_rlsoda", PACKAGE = "virosolver")
}
