#' Example SEIR model parameter table
#'
#' Example of the parameter table used for to simulate from an SEIR model with \code{virosolver}. 
#' This data frame is used to control everything related to the model parameters, including values, names, fixed/estimates, and uniform prior bounds.
#' For this parameter table, we have entries for the SEIR model parameters (R0, infectious period, incubation period, seed time and seed size), as well as
#' entries for the viral kinetics/Ct model (peak Ct value, waning duration, time to peak, variation about the model Ct curve etc).
#' @docType data
#' @usage data(example_seir_partab)
#' @format A data frame with 20 rows and 8 variables:
#' \describe{
#'     \item{values}{numeric values of the parameters}
#'     \item{names}{string names of the model parameters}
#'     \item{fixed}{binary values indicating if the parameter should be fixed (1) or estimated (0) during the MCMC procedure}
#'     \item{lower_bound}{lower numeric bound for the parameter during fitting (lower uniform prior bound)}
#'     \item{upper_bound}{upper numeric bound for the parameter during fitting (upper uniform prior bound)}
#'     \item{steps}{value between 0 and 1, giving the initial step size in the MCMC proposals. Note that these are adapted automatically}
#'     \item{lower_start}{can be used to set the lower allowable random starting value for the MCMC}
#'     \item{upper_start}{can be used to set the upper allowable random starting value for the MCMC}
#' }
#' @family example_data
"example_seir_partab"

#' Example Gaussian Process prior parameter table
#'
#' Example of the parameter table used for to fit the Gaussian Process prior model with \code{virosolver}. 
#' This data frame is used to control everything related to the model parameters, including values, names, fixed/estimates, and uniform prior bounds.
#' For this parameter table, we have entries for the GP prior parameters (the latent states for each day "prob", the hyperparameters for the exponential kernel function, rho and nu), as well as
#' entries for the viral kinetics/Ct model (peak Ct value, waning duration, time to peak, variation about the model Ct curve etc).
#' @docType data
#' @usage data(example_gp_partab)
#' @format A data frame with 384 rows and 8 variables:
#' \describe{
#'     \item{values}{numeric values of the parameters}
#'     \item{names}{string names of the model parameters}
#'     \item{fixed}{binary values indicating if the parameter should be fixed (1) or estimated (0) during the MCMC procedure}
#'     \item{lower_bound}{lower numeric bound for the parameter during fitting (lower uniform prior bound)}
#'     \item{upper_bound}{upper numeric bound for the parameter during fitting (upper uniform prior bound)}
#'     \item{steps}{value between 0 and 1, giving the initial step size in the MCMC proposals. Note that these are adapted automatically}
#'     \item{lower_start}{can be used to set the lower allowable random starting value for the MCMC}
#'     \item{upper_start}{can be used to set the upper allowable random starting value for the MCMC}
#' }
#' @family example_data
"example_gp_partab"


#' Example Ct value data
#'
#' Observed Ct values simulated from an SEIR model, with 2000 individuals sampled at random every 14 days between days 55 and 153 of the simulation. Each entry is one anonymised individual sample.
#' @docType data
#' @usage data(example_ct_data)
#' @format A data frame with 15984 rows and 2 variables:
#' \describe{
#'     \item{t}{Time in days that the sample was taken, relative to the start time of the simulation (simulation start time on t=0)}
#'     \item{ct}{Ct value of the sample}
#' }
#' @family example_data
"example_ct_data"

#' True simulated SEIR incidence
#' 
#' The true SEIR incidence curve used in the case study simulation
#' @docType data
#' @usage data(example_seir_incidence)
#' @format A tible with 251 rows and 2 columns
#' \describe{
#'     \item{t}{Time in days}
#'     \item{prob_infection}{Per capita incidence}
#' }
#' @family example_data
"example_seir_incidence"

#' Simulated SEIR incidence for vignettes
#' 
#' This is the raw incidence from the simulate_seir_process function. It is created in vignette 1 and used in vignette 2. 
#' @docType data
#' @usage data(vignette_data1)
#' @format Numeric vector with length 151
#' @family example_data
"vignette_data1"

#' Simulated Ct values for vignettes
#' 
#' These are the simulated Ct values which are created in vignette 1 and used in vignette 2. 
#' @docType data
#' @usage data(vignette_data2)
#' @format A data frame with 22237 rows and and 4 columns
#' \describe{
#'     \item{sampled_time}{Time of sample collection in days}
#'     \item{Ct}{Ct values}
#'     \item{date}{calendar date}
#'     \item{Ct_round}{rounded Ct values}
#' }
#' @family example_data
"vignette_data2"
