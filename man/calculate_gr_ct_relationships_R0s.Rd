% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maths_functions.R
\name{calculate_gr_ct_relationships_R0s}
\alias{calculate_gr_ct_relationships_R0s}
\title{Calculate GR CT relationship for multiple R0s}
\usage{
calculate_gr_ct_relationships_R0s(
  model_pars,
  viral_load_pars,
  times,
  lastday = 35,
  cts = seq(0, 40, by = 0.1),
  R0s,
  symptom_surveillance = FALSE
)
}
\description{
Calculate GR CT relationship for multiple R0s
}
